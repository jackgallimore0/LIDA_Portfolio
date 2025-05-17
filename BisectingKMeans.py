import random
import math
import numpy as np
import matplotlib.pyplot as plt
import os

random.seed(42)

# Loads dataset from a file into a list of tuples (index, vector)
# Input: filename (string, default "dataset")
# Output: list of tuples [(index, vector), ...]
def load_dataset(filename="dataset"):
    try:
        with open(filename, "r") as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return []
    except IOError:
        print(f"Error: Could not read file '{filename}'.")
        return []

    X = []
    for i, line in enumerate(lines):
        words = line.strip().split()
        try:
            tup = (i, [float(word) for word in words[1:]])
            X.append(tup)
        except ValueError:
            print(f"Warning: Skipping corrupted line {i + 1}.")

    return X

# Chooses initial centroids for k-means++ algorithm
# Inputs: 
# - X: dataset (list of tuples)
# - k: number of clusters (int)
# Output: list of initial centroids [(index, vector), ...]
def initialSelection(X, k):
    Y = [random.choice(X)] # Choose the first centroid randomly from X
    
    for _ in range(1, k):
        distances = []
        # For each data point, compute the squared distance to the nearest chosen centroid.
        for x in X:
            # point is (index, vector)
            index, vector= x
            # Calculate squared distances to all centroids in Y
            sq_dists = [sum((x - y) ** 2 for x, y in zip(vector, center[1])) for center in Y]
            #print(sq_dists)
            distances.append(min(sq_dists))
            #print(f"minimum squared distance {min(sq_dists)}")
            
        total_distance = sum(distances)
        r = random.uniform(0, total_distance)
        cumulative = 0
        
        # Choose the new centroid weighted by the squared distances
        for point, d in zip(X, distances):
            cumulative += d
            if cumulative >= r:
                Y.append(point)
                break
                
    return Y

# Calculates Euclidean distance between two vectors
# Input: vectors a, b
# Output: float distance
def eucl_distance(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))

# Assigns data points to the nearest centroid
# Inputs: 
# - X: dataset (list of tuples)
# - k: number of clusters
# - Y: centroids (list of tuples)
# Output: clusters (list of lists of tuples)
def assignClusterIds(X, k, Y):
    #Split X into a list of lists
    #Each list is a list of clusters [[(point,vector),(point,vector),....],[(point,vector),....],.....]
    clusters = [[] for _ in range(k)]

    for x in X:
        index, vector = x
        distances = [eucl_distance(vector, centroid[1]) for centroid in Y]
        min_index = distances.index(min(distances))
        clusters[min_index].append((index, vector))
    return clusters

# Computes new centroid positions by averaging points in each cluster
# Input: clusters (list of lists of tuples)
# Output: new centroids (list of tuples)
def computeClusterRepresentatives(clusters):
    #Creates new Y same format as Y
    #Uses index -1 as they are not original points
    new_Y = []
    for cluster in clusters:
        vector_list = []
        for x in cluster:
            _, vector = x
            vector_list.append(vector)
        new_mean = list(np.mean(vector_list, axis=0))
        new_Y.append((-1, new_mean))
    return new_Y

# Computes silhouette score to evaluate clustering quality
# Input: clusters (list of clusters)
# Output: mean silhouette coefficient (float)
def computeSilhouette(clusters):
    sihlohettes = [] #List of tuples (index, score)
    
    for i, cluster_i in enumerate(clusters):
        a = 0
        b = math.inf #inf b as we need to find b
        #Comapares each cluster to each cluster including itself
        for j, cluster_j in enumerate(clusters):
            total_distance = 0
            count = 0
            for _, vec_i in cluster_i:
                for _, vec_j in cluster_j:
                    #Sum of euclidian distance
                    dist = eucl_distance(vec_i, vec_j)
                    total_distance += dist
                    count += 1
                    
            average_distance = total_distance / count
            #print(f"Average distance between cluster {i} and cluster {j}: {average_distance:.4f}")
            #if i = j it is a
            #if i =/ j we compare it to currrent b to see if it is lower
            if i == j:
                a = average_distance
            if i != j:
                if average_distance < b:
                    b = average_distance

        #calculate sihloette, if only 1 cluster b remains asi inf thus sihloette = 0   
        if b == math.inf:
            shilohette = 0 
        else:
            shilohette = (b - a) / max(a, b)
                    
        #print(f"shiloette for: {i} is: {shilohette}")
        sihlohettes.append((i, shilohette))
        
    mean_sihlohette = sum(x[1] for x in sihlohettes) / len(sihlohettes)
    return (mean_sihlohette)

# Plots silhouette scores against the number of clusters
# Input: silhouette_scores (list of floats)
# Output: None (plot saved to file)
def plotSihlohette (silhouette_scores):
    #Plot the shiloette scores
    plt.figure(figsize=(8, 5))
    k_values = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    plt.bar(k_values, silhouette_scores)
    plt.title("Silhouette Coefficient vs. Number of Clusters (k)")
    plt.xlabel("Number of Clusters (k)")
    plt.ylabel("Silhouette Coefficient")
    plt.xticks(k_values)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.savefig("BisectingKMenas Plot", bbox_inches='tight')
    plt.close()

# Runs k-means with k=2, used for bisecting k-means
# Input: X (dataset)
# Output: clusters and centroids
def Kmeans2 (X):
    #Runs the process of splitting a cluster into 2 clusters
    Y = initialSelection(X, 2)
    clusters = assignClusterIds(X, 2, Y)
    new_Y = computeClusterRepresentatives(clusters)
    return(clusters, new_Y)

# Calculates total distances within clusters to find the one to split
# Inputs:
# - clusters (list of clusters)
# - Y (list of centroids)
# Output: list of summed distances for each cluster
def SumDistSqu (clusters, Y):
    #returns list of squared distances of clusters s we can choose ma to split
    sumDistlist = []
    for i in range(len(clusters)):
        cluster = clusters[i]
        _, y_vector = Y[i]
        sumDist = 0
        
        for point in cluster:
            _, x_vector = point
            #print(f"x vector:{x_vector}")
            #print(f"y_vector:{y_vector}")
            sumDist += eucl_distance(x_vector,y_vector)
        
        sumDistlist.append(sumDist)
    return sumDistlist

#Main
#You don’t need to run bisecting for every value of s, since the process already follows the path that gives the correct silhouette coefficient at each level.
X = load_dataset()
silhouette_scores = [0]

clusters, Y = Kmeans2(X) #Splits in 2
silhouette_scores.append(computeSilhouette(clusters))

for i in range (2,9):
    sumdistlsit = SumDistSqu(clusters, Y)
    max_index = sumdistlsit.index(max(sumdistlsit))
    expand_cluster = clusters[max_index]
    clusters.pop(max_index)
    Y.pop(max_index)
    expanded_cluster, new_Y  = Kmeans2(expand_cluster)
    clusters.extend(expanded_cluster)
    Y.extend(new_Y)
    silhouette_scores.append(computeSilhouette(clusters))

plotSihlohette(silhouette_scores)
