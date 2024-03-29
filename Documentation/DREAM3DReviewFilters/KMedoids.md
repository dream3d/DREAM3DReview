# K Medoids #

## Group (Subgroup) ##

DREAM3D Review (Clustering)

## Description ##

This **Filter** applies the k medoids algorithm to an **Attribute Array**.  K medoids is a _clustering algorithm_ that assigns to each point of the **Attribute Array** a _cluster Id_.  The user must specify the number of clusters in which to partition the array.  Specifically, a k medoids partitioning is such that each point in the data set is associated with the cluster that minimizes the sum of the pair-wise distances between the data points and their associated cluster centers (medoids).  This approach is analogous to [k means](@ref kmeans), but uses actual data points (the medoids) as the cluster exemplars instead of the means.  Medoids in this context refer to the data point in each cluster that is most like all other data points, i.e., that data point whose average distance to all other data points in the cluster is smallest.  Unlike [k means](@ref kmeans), since pair-wise distances are minimized instead of variance, any arbirtary concept of "distance" may be used; this **Filter** allows for the selection of a variety of distance metrics.    

This **Filter** uses the _Voronoi iteration_ algorithm to produce the clustering.  The algorithm is iterative and proceeds as follows:

1. Choose k points at random to serve as the initial cluster medoids
2. Associate each point to the closest medoid
3. Until convergence, repeat the following steps:
  * For each cluster, change the medoid to the point in that cluster that minimizes the sum of distances between that point and all other points in the cluster
  * Reassign each point to the closest medoid

Convergence is defined as when the medoids no longer change position.  Since the algorithm is iterative, it only serves as an approximation, and may result in different classifications on each execution with the same input data.  The user may opt to use a mask to ignore certain points; where the mask is _false_, the points will be placed in cluster 0.
    
A clustering algorithm can be considered a kind of segmentation; this implementation of k medoids does not rely on the **Geometry** on which the data lie, only the _topology_ of the space that the array itself forms.  Therefore, this **Filter** has the effect of creating either **Features** or **Ensembles** depending on the kind of array passed to it for clustering.  If an **Element** array (e.g., voxel-level **Cell** data) is passed to the **Filter**, then **Features** are created (in the previous example, a **Cell Feature Attribute Matrix** will be created).  If a **Feature** array is passed to the **Filter**, then an **Ensemble Attribute Matrix** is created.  The following table shows what type of **Attribute Matrix** is created based on what sort of array is used for clustering:

| Attribute Matrix Source             | Attribute Matrix Created |
|------------------|--------------------|
| Generic | Generic |
| Vertex | Vertex Feature | 
| Edge | Edge Feature |
| Face | Face Feature | 
| Cell | Cell Feature| 
| Vertex Feature | Vertex Ensemble |
| Edge Feature | Edge Ensemble |
| Face Feature | Face Ensemble |
| Cell Feature | Cell Ensemble|
| Vertex Ensemble | Vertex Ensemble |
| Edge Ensemble | Edge Ensemble |
| Face Ensemble | Face Ensemble |
| Cell Ensemble | Cell Ensemble|

This **Filter** will store the medoids for the final clusters within the created **Attribute Matrix**.

## Parameters ##

| Name | Type | Description |
|------|------|-------------|
| Number of Clusters | int32_t | The number of clusters in which to partition the array |
| Distance Metric | Enumeration | The metric used to determine the distances between points |
| Use Mask | bool | Whether to use a boolean mask array to ignore certain points flagged as _false_ from the algorithm |
| Use Random Seed | bool | Use a user defined random seed value |
| Random Seed Value | uint64 | The random seed to use |

## Required Geometry ###

None

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| Any **Attribute Array** | None | Any| Any | The **Attribute Array** to cluster |
| **Attrubute Array** | Mask | bool | (1) | Specifies if the point is to be counted in the algorithm, if _Use Mask_ is checked |

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Attribute Matrix** | ClusterData | Feature/Ensemble | N/A | The **Attribute Matrix** in which to store information associated with the created clusters |
| **Attribute Array** | ClusterIds | int32_t | (1) | Specifies to which cluster each point belongs |
| **Attribute Array** | ClusterMeans | double | (1) | The means of the final clusters |


## References ## 

[1] A simple and fast algorithm for K-medoids clustering, H.S. Park and C.H. Jun, Expert Systems with Applications, vol. 28 (2), pp. 3336-3341, 2009.

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users
