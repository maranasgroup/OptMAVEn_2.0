import os
from collections import OrderedDict

import numpy as np
from sklearn.decomposition import PCA

import standards


def kmeans(data, k, tolerance=standards.DefaultKmeansTolerance, max_iterations=standards.DefaultKmeansMaxIterations):
    if k > len(data):
        raise ValueError("Cannot set k to {}, greater than number of data ({}).".format(k, len(data)))
    centroids = random_centroids(data, k)  
    movement = None
    iteration = 1
    while iteration <= max_iterations and (movement is None or movement > tolerance):
        clusters = update_clusters(data, centroids)
        movement = update_centroids(clusters, centroids)
        iteration += 1
    #FIXME did I correctly implement this error calculation?
    max_mse = np.max([np.mean((np.array([to_array(vector) for vector in cluster]) - centroid)**2) for cluster, centroid in zip(clusters, centroids)])
    return clusters, max_mse
   

def update_clusters(data, centroids):
    k = len(centroids)
    clusters = [list() for i in range(k)]
    for item in data:
        nearest_cluster_index = np.argmin([np.linalg.norm(to_array(item) - centroid) for centroid in centroids])
        clusters[nearest_cluster_index].append(item)
    # If any cluster is empty then assign one point
    # from data set randomly so as to not have empty
    # clusters and 0 means.
    for i in range(k):
        cluster_i = clusters[i]
        if len(cluster_i) == 0:
            for j in np.random.permutation(range(k)):
                cluster_j = clusters[j]
                n_j = len(cluster_j)
                if n_j > 1:
                    cluster_i.append(cluster_j.pop(np.random.randint(n_j)))
    return clusters


def update_centroids(clusters, centroids):
    if len(clusters) != len(centroids):
        raise ValueError("Clusters and Centroids have different lengths.")
    square_distance = 0.0
    for index, cluster in enumerate(clusters):
        displacement = np.mean([to_array(vector) for vector in cluster], axis=0) - centroids[index]
        square_distance += np.dot(displacement, displacement)
        centroids[index] += displacement
    return np.sqrt(square_distance)


def random_centroids(data, k):
    indexes = np.random.choice(range(len(data)), size=k, replace=False)
    centroids = [to_array(data[index]) for index in indexes]
    return centroids


def to_array(iterable):
    if isinstance(iterable, np.ndarray):
        return iterable
    else:
        return np.array(list(iterable))


def normalize(vector, means, stds):
    """ Return a normalized vector by subtracting the means and dividing by the std. devs. """
    return np.nan_to_num((to_array(vector) - means) / stds)


def calculate_means_stds(data):
    """ Calculate the mean and standard deviation of each feature in the data. """
    d_array = np.array([to_array(item) for item in data])
    return np.mean(d_array, axis=0), np.std(d_array, axis=0)


def pca_run(data):
    """ Perform PCA. """
    means, stds = calculate_means_stds(data)
    d_array = normalize(np.array([to_array(item) for item in data]), means, stds)
    pca = PCA(n_components=3)
    d_3d = pca.fit_transform(d_array)
    return [tuple(x) for x in d_3d]
    

def nopca_run(data):
    """ Do not perform PCA; nopca_run exists to minic the interface of the pca_run function. """
    means, stds = calculate_means_stds(data)
    d_array = normalize(np.array([to_array(item) for item in data]), means, stds)
    return [tuple(to_array(item)) for item in d_array]


def optimal_kmeans(data, threshold=standards.DefaultKmeansOptKThreshold, tolerance=standards.DefaultKmeansTolerance, max_iterations=standards.DefaultKmeansMaxIterations):
    data_pc = pca_run(data)
    pc_map = {pc: datum for datum, pc in zip(data, data_pc)}  # map the PC-transformed data back to the original data objects
    k = 1 
    mses = OrderedDict()
    # Keep incrementing k until an optimal k level is found.
    while k <= len(data_pc) and (k == 1 or min(mses.values()) / mses[1] > threshold):
        print("\nK =", k)
        print(data_pc)
        clusters, mse = kmeans(data_pc, k, tolerance, max_iterations)
        if k == 1 or mse < min(mses.values()):
            clusters_opt = clusters
        mses[k] = mse
        print("MSE", mse)
        print(min(mses.values()) / mses[1], threshold)
        k += 1
    # perform the mapping
    clusters_opt = [[pc_map[datum_pc] for datum_pc in cluster] for cluster in clusters_opt]
    return clusters_opt, mses
