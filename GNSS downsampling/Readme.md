# Dsample.m (Downsampling of GNSS or seismic stations)

In any case, where downsampling of stations is necessary for further calculation.  
Assures the distance between stations is under a self-defined distance.  
Remove one station at a time.

1. Determine a distance threshold for which stations should be spaced
2. Calculate the number of stations that are within the distance for each station
3. Find the stations that have the most neighboring stations
4. Calculate the total distance between these stations
5. Remove the station that has the minimum total distance (Closest to other stations)

##### DistMatrix.m is a sub-routine in Dsample.m. Must be included in path.

---
##### Input variable:
   * InMat: matrix. The first two columns have to be local coordinates. 
   * Dist: number. Distance threshold between stations.
##### Output variable:
   * OutMat: matrix. Downsampled data matrix.
   * OmitInd: number. The index of downsampled stations.
---
### Example:
```MatLab
% Determine distance threshold (Meters)
Dist = 5000;

% Start program Dsample.m
[Downsampled,NInd] = Dsample(ToBeDownsampled,Dist);
```
After 272 iterations (576 stations downsampled to 304 under distance threshold of 5000 meters):
![Example](https://github.com/LiChiehLin/3D_decomposition/blob/43852c0bc77ba0cb97f794eeb8fc51ee38bf16ce/Figures/Dsample_Example.png)
