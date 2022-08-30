##### Input variables:
- InGrd: 
   * cell variable. Each column in the cell should be a image matrix read by *'grdread2.m'*
- Azimuth: 
   * Satellite inclination degree
   * Ascending and Descending should be the same number
   * 12° for Sentinel-1 
   * 8° for ALOS-2
- LookAngle:
   * Satellite looking angle
   * 32.9° for Sentinel-1 subswath 1 
   * 38.3° for Sentinel-1 subswath 2 
   * 43.1° for Sentinel-1 subswath 3
- DispType:
   * cell variable. Each input is a string either *'LOS'* or *'Azi'*
- Orbit: 
   * cell variable. Each input is a string either *'Asc'* or *'Des'*
- Num:
   * number. How many images are under converting


### Note that:
Each input variable should align to the image attributes  
