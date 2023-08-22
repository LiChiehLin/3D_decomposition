# HFsep.m (Hanging and footwall separation)

For some post-processing where one desires to separate the hanging wall and footwall into two matrix blocks. (i.e. Denoise)  
This code takes on the distinct different displacement pattern across two fault blocks.  
Works more satisfactorily on linear fault trace.  
You can always use GIS software to mask out the blocks, this is just an alternative.  

##### findlocalmin.m and findlocalmax.m are sub-routines in HFsep.m. Must be included in path.

---
##### Input variable:
   * InGrd: matrix. Matrix that needs to be separated.
   * Strike: number. The striking direction of the fault. Direction respect to north.  
     **Support input degree ranging from 0 to 360 (decimal number is accepted)**.  
     - North-South fault: input 0 or 180 degrees.  
     - Northeast-Southwest fault: input 45 or 225 degrees.  
     - East-West fault: input 90 or 270 degrees.
   * PtoN: number. Matrix value going from positive to negative or not.
     - 1: Yes; 0: No
  * Overlap: number: How much overlapping is there between the two blocks.  
##### Output variable:
   * SideA: matrix. One side of the matrix.
   * SideB: matrix. The other side of the matrix.
   * SideAIndex: vector. The index of SideA matrix.
   * SideBIndex: vector. The index of SideB matrix.

Use `MatComb.m` to combine the separated matrix into one matrix.

---  
### Example:
```MatLab
% Load displacement matrix
[x,y,z] = grdread2('Sentinel_RangeOffset.grd');
% Determine parameters for HFsep.m
Str = 93.1;
Overlap = 1;
PtoN = 0;

% Start program Dsample.m
[Hang,Foot,HangInd,FootInd] = HFsep(z,Str,PtoN,Overlap);
```
