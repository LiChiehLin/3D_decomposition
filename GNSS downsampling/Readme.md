# Dsample.m (Downsampling of GNSS or seismic stations)

In any case, where downsampling of stations is necessary for further calculation. 
Assures the distance between stations is under a self-defined distance.

1. Determine a distance threshold for which stations should be spaced
2. Calculate the number of stations that are within the distance for each station
3. Calculate the Cumulative Distribution Function (CDF) of all `σ`
4. Determine `Threshold γ` based on 95% of the CDF
5. Remove pixels that the corresponding `σ` is larger than `γ`

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
% Read data and determine window size
[x,y,ToBeRemoved] = grdread2('Displacement.grd');
ws = 3;

% Start program Deno.m
[Removed,NInd] = Deno(ToBeRemoved,ws);
```
