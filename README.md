# 3D_decomposition
> This repository contains a few MatLab codes for converting SAR-derived displacement into true 3D displacement (EW, NS and UD directions)  
> The code for reading netCDF grds is from Kelsey Jordahl (2022). grdread2 (https://www.mathworks.com/matlabcentral/fileexchange/25683-grdread2), MATLAB Central File Exchange. Retrieved August 30, 2022.  

- Read grd: grdread2.m  
- 3D decomposition program: InSAR3Ddisp.m  
- Remove noisy pixels: Deno.m  
- Convert to netCDF format: Convert2Grd.m  
- Get image boundaries and increments: BoundBox.m  
- Regrid the grds to the same size: Regrid.m  
- Sptial fitering: filtsp.m  

