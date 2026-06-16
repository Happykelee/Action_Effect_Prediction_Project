"""
Paired Spatio-Temporal Cluster-Based Permutation Testing
========================================================
Custom implementation for M/EEG data based on Mike X Cohen's methodology (Cohen, 2017 p.233-244).
The multiple comparison correction is based on cluster size and mass under the null hypothesis.

Dependencies: random, numpy, scipy, joblib
"""

import random
import numpy as np
import pandas as pd
from scipy import stats
from joblib import Parallel, delayed

class UnionFind:
    """
    Union-Find (Disjoint Set) data structure.
    Used to efficiently track and merge connected components in a graph. Here, it is applied 
    to rapidly identify spatiotemporally connected significant data points (i.e., clusters).
    """
    def __init__(self): 
        '''Dictionary to store the parent of each node.'''
        self.parent = {}
        
    def add(self, x):
        """Adds a new node, defaulting its parent to itself."""
        if x not in self.parent: 
            self.parent[x] = x
            
    def find(self, x):
        """
        Finds the root node of the set containing the node.
        Incorporates Path Compression to accelerate future lookups.
        """
        self.add(x)
        if self.parent[x] != x: 
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]
        
    def union(self, x, y):
        """Merges the sets containing two nodes."""
        self.add(x); self.add(y)
        rx = self.find(x); ry = self.find(y)
        if rx != ry: 
            self.parent[ry] = rx

def find_connected_regions(binary_matrix, topology_map):
    """
    Identifies connected regions (clusters) within a binary threshold matrix.
    
    Parameters:
    binary_matrix: ndarray, a boolean matrix of shape (n_channels, n_times). 
                   True indicates the statistic exceeded the uncorrected threshold.
    topology_map: list of arrays, the topological adjacency map of electrodes. 
                  topology_map[i] contains indices of electrodes spatially adjacent to electrode i.
    
    Returns:
    clusters: list of lists, all identified clusters. 
              Each cluster is a list of (channel_idx, time_idx) tuples.
    """
    # Extract electrode and time indices of all significant points
    e_idx, t_idx = np.where(binary_matrix)
    points = list(zip(e_idx, t_idx))
    if not points: 
        return []
    
    point_set = set(points)
    max_t = binary_matrix.shape[1]
    uf = UnionFind()
    
    for pt in points:
        e, t = pt
        # 1. Check for spatial adjacency (Spatial neighbors across electrodes)
        for n_e in topology_map[e]:
            if (n_e, t) in point_set: 
                uf.union(pt, (n_e, t))
                
        # 2. Check for temporal adjacency (Temporal neighbors, +/- 1 time point)
        for dt in (-1, 1):
            if 0 <= t + dt < max_t and (e, t + dt) in point_set: 
                uf.union(pt, (e, t + dt))
            
    # Group points belonging to the same root node (i.e., the same cluster)
    clusters = {}
    for pt in points:
        root = uf.find(pt)
        if root not in clusters: 
            clusters[root] = []
        clusters[root].append(pt)
        
    return list(clusters.values())

def execute_single_permutation(diff_tensor, n_subs):
    """
    Executes a single permutation operation.
    Under the null hypothesis, the signs of the pre-post differences are randomly flipped.
    
    Parameters:
    diff_tensor: ndarray, difference matrix of shape (n_subs, n_times, n_channels).
    n_subs: int, number of subjects.
    
    Returns:
    t_stat: ndarray, pseudo T-statistic distribution map of shape (n_times, n_channels).
    """
    # Randomly generate a sign-flip vector (half 1, half -1) and shuffle it
    sign_flips = np.r_[np.ones(n_subs // 2), -np.ones(n_subs - n_subs // 2)]
    random.shuffle(sign_flips)
    
    # Apply the flip vector to the difference tensor
    permuted_diff = diff_tensor * sign_flips[:, np.newaxis, np.newaxis]
    
    # Recalculate the one-sample T-test (against 0) for the permuted data
    return stats.ttest_1samp(permuted_diff, 0).statistic

def extract_max_null_clusters(perm_t_map, cluster_forming_thresh, topology_map):
    """
    Extracts the maximum cluster features (size and mass) from a single permuted T-statistic map.
    This builds the null distribution required to control the Family-Wise Error Rate (FWER).
    
    Parameters:
    perm_t_map: ndarray, the permuted T-statistic map.
    cluster_forming_thresh: float, initial T-value threshold for forming clusters.
    topology_map: list, spatial adjacency map of electrodes.
    
    Returns:
    null_metrics: ndarray, an array of 4 elements: 
                  [Max Positive Size, Max Negative Size, Max Positive Mass, Max Negative Mass]
    """
    null_metrics = np.zeros(4) 
    
    # Process positive significant (> thresh) and negative significant (< -thresh) masks
    for i, mask in enumerate([perm_t_map > cluster_forming_thresh, perm_t_map < -cluster_forming_thresh]):
        found_clus = find_connected_regions(mask.T, topology_map)
        
        if found_clus:
            # Calculate the sizes (number of data points) of all found clusters
            sizes = [len(c) for c in found_clus]
            max_idx = np.argmax(sizes) # Index of the largest cluster
            
            # Record the size of the maximum cluster
            null_metrics[i] = sizes[max_idx]
            # Record the mass/sumstats of the maximum cluster (sum of all T-values within it)
            null_metrics[i + 2] = np.sum([perm_t_map[j, k] for k, j in found_clus[max_idx]])
            
    return null_metrics

def run_cohen_spatiotemporal_permutation(diff_tensor, n_permutations=10000, alpha_forming=0.02,
                                         topology_map=None, random_seed=0):
    """
    Main function to run the spatiotemporal permutation test with cluster-based correction based on Mike X Cohen's methodology (Cohen, 2017 p.233-244).
    
    Parameters:
    diff_tensor: ndarray, difference matrix of shape (n_subs, n_times, n_channels).
    n_permutations: int, number of permutations (e.g., 10000).
    alpha_forming: float, initial alpha threshold for forming clusters (uncorrected p-value, e.g., 0.02).
    topology_map: list, adjacency map of the electrodes.
    random_seed: int, base seed for parallel processing to ensure reproducibility.
    
    Returns:
    obs_t_map: ndarray, T-statistic map of the actual observed data.
    all_clusters: list, coordinates of all clusters found in the actual data.
    summary: DataFrame, contains cluster sizes, masses, and their non-parametric p-values.
    null_distributions: ndarray, the generated null distributions of shape (n_permutations, 4).
    """
    random.seed(a=random_seed)
    n_subs = diff_tensor.shape[0]
    # Calculate the cluster-forming T-value threshold (two-tailed)
    thresh_t = stats.t.ppf(1 - alpha_forming / 2, n_subs - 1)
    
    # 1. Compute observed statistics
    obs_t_map = stats.ttest_1samp(diff_tensor, 0).statistic
    real_pos_clusters = find_connected_regions((obs_t_map > thresh_t).T, topology_map)
    real_neg_clusters = find_connected_regions((obs_t_map < -thresh_t).T, topology_map)
    all_clusters = real_pos_clusters + real_neg_clusters
    
    # 2. Parallel permutation execution 
    perm_maps = Parallel(n_jobs=-1)(
        delayed(execute_single_permutation)(diff_tensor, n_subs) for i in range(n_permutations)
    )
    
    # 3. Build null distributions
    null_distributions = Parallel(n_jobs=-1)(
        delayed(extract_max_null_clusters)(pm, thresh_t, topology_map) for pm in perm_maps
    )
    null_distributions = np.stack(null_distributions)
    
    # 4. Calculate cluster sizes, masses, and their non-parametric p-values
    summary = {'size': [], 'sumstats': [], 'pval_size': [], 'pval_sumstats': []}
    for c_idx, clusters in enumerate([real_pos_clusters, real_neg_clusters]):
        for clu in clusters:
            sz = len(clu)
            # Calculate the mass (sum of T-values) for the current real cluster
            mass = np.sum([obs_t_map[j, k] for k, j in clu])
            summary['size'].append(sz)
            summary['sumstats'].append(mass)
            # P-value based on cluster size: proportion of null max sizes > real size
            summary['pval_size'].append(
                (np.sum(null_distributions[:, c_idx] > sz) + 1) / (n_permutations + 1))
            # P-value based on cluster mass: proportion of null max masses > real mass
            summary['pval_sumstats'].append(
                (np.sum(np.abs(null_distributions[:, c_idx + 2]) > np.abs(mass)) + 1) / (n_permutations + 1))
            
    return obs_t_map, all_clusters, pd.DataFrame(summary), null_distributions