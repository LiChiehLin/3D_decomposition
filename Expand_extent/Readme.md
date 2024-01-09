# Expand_extent.m 

Expand the spatial extent of the input grds. This is convenient when performing matrix subtraction or addition. The expanded grids will be filled with *NaN*.

### Note that the x and y increments should be the same to prevent grid value mis-registration. Therefore, use with cautious.

---
##### Input variable:
   * InMat: cell. Containing the Lon, Lat and value of each input grd (format is the one `grdread2.m` records)  
     * col1: Lon vector
     * col2: Lat vector
     * col3: z-value
##### Output variable:
   * OutMat: cell. Containing the expaned Lon, Lat and grds.
     * row1: Expanded Lon (Can be used in `grdwrite2.m`)
     * row2: Expanded Lat (Can be used in `grdwrite2.m`)
     * row3~end: Expanded input grds 
---
### Example:
```MatLab
% Example: Two grds need to be expanded
[x1,y1,z1] = grdread2('grd1.grd');
[x2,y2,z2] = grdread2('grd2.grd');
InMat{1,1} = x1; InMat{1,2} = y1; InMat{1,3} = z1;
InMat{2,2} = x2; InMat{2,2} = y2; InMat{2,3} = z2;

% Start program Extract_Timespan.m
OutMat = Expand_extent(InMat);
```
#### Output as grd
```MatLab
grdwrite2(OutMat{1},OutMat{2},OutMat{3});
grdwrite2(OutMat{1},OutMat{2},OutMat{4});
```
