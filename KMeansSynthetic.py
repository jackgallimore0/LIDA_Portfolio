import random
import math
import numpy as np
import matplotlib.pyplot as plt
import os

random.seed(42)

# Loads data from a file into a list of tuples.
# Input: filename (string) - path to the data file.
# Output: list of tuples (index, vector).
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

# Selects initial random centroids from data points.
# Inputs:
# - X: dataset (list of tuples)
# - k: number of clusters (int)
# Output: list of tuples representing centroids (index, vector).
def initialSelection (X, k):
    #Takes k random data points as representatives list of (index, data)
    random_locations = random.sample(range(len(X)), k)
    random.seed(42)
    Y = [X[i] for i in random_locations]
    return(Y)

# Calculates Euclidean distance between two vectors.
# Inputs: two vectors a, b
# Output: float representing the distance.
def eucl_distance(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))

# Assigns data points to clusters based on nearest centroid.
# Inputs:
# - X: dataset (list of tuples)
# - k: number of clusters (int)
# - Y: centroids (list of tuples)
# Output: list of clusters, each cluster is a list of tuples.
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

# Computes new centroids by averaging cluster points.
# Input: clusters (list of clusters)
# Output: list of tuples representing new centroids.
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

# Computes the silhouette score to evaluate clustering quality.
# Input: clusters (list of clusters)
# Output: float (average silhouette score).
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

# Plots silhouette scores vs. number of clusters and saves to file.
# Input: silhouette_scores (list of floats)
# Output: None (plot saved as file).
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

    plt.savefig("KMenasSynthetic Plot", bbox_inches='tight')
    plt.close()

#Main
silhouette_scores = []
X = [(i, [random.uniform(-1, 1) for _ in range(300)]) for i in range(327)]

for k in range(1, 10):
    Y = initialSelection(X, k)
    clusters = assignClusterIds(X, k, Y)
    new_Y = computeClusterRepresentatives(clusters)
    s = computeSilhouette(clusters)
    silhouette_scores.append(s)
    
plotSihlohette(silhouette_scores)

