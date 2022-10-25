# nearneighbor.m (Interpolation)

Fill NaN pixels with values from neighboring pixels  

1. How many neighbors are searched is defined by `Radius`  
  If the NaN pixel is `(i,j)` and `Radius` is set to 3  
  Then the neighbors are given by `(i-Radius:i+Radius,j-Radius:j+Raidus)`
  
2. Weighting of neighboring pixels is defined as (Radius is 3)  
  Say neighbors in `Raidus 1` is `a`, `Raidus 2` is b and `Radius 3` is `c`  
  The interpolated-value will be determined as `mean(a+(a+b)+(a+b+c))`  
  So the weighting will be 3:2:1


---
##### Input variable:
   * InMat: matrix. Read from *grdread2.m*
   * Radius: number. The radius of neighbors to be searched
   * Weight: number. Whether or not give weighting to neighboring values (1:Yes, 0:No)
##### Output variable:
   * OutMat: matrix. Interpolated image matrix
---
### Example:
```MatLab
% Read data and determine window size
[x,y,ToBeInterpolated] = grdread2('Displacement.grd');
Radius = 3;
Weight = 1

% Start program nearneighbor.m
Interpolated = nearneighbor(ToBeInterpolated,Radius,Weight);
```

