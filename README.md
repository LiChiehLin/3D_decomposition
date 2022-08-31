# 3D_decomposition
This repository contains a few MatLab codes for converting SAR-derived displacement into true 3D displacement (EW, NS and UD directions) and some useful tools to handle displacement data  

These codes were developed under MatLab 2015b  
  
Put every code in one directory and add path to your *main* script.  
```MatLab
addpath(genpath('DIRECTORY OF MATLAB CODES'))
```

---
Use grdread2.m to read netCDF grds
> The code for reading netCDF grds is from Kelsey Jordahl (2022). grdread2 (https://www.mathworks.com/matlabcentral/fileexchange/25683-grdread2), MATLAB Central File Exchange. Retrieved August 30, 2022.  

### Below lists the codes in each directories
- 3D decomposition: 
   * InSAR3Ddisp.m  
   * Convert2Grd.m
- Denoise: Deno.m  
- Sptial fitering: filtsp.m  
