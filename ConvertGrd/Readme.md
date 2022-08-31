# Convert2Grd

* Inv: cell variable
   - (Lon Lat E N U) or (Lon Lat Z) matices with same sizes 
* Inc: vector variable. X increment and Y increment
   - product from *BoundBox.m* Or you can make by yourself
* Step: string variable
   - The name of the converted grd



### Example (3D displacement field):
All matrices (`Lon` `Lat` `E` `N` `U`) should be the same size (M-by-N):

```MatLab
[x,y,z] = grdread2('Sentinel_Asc_LOS.grd');
[x2,y2,z2] = grdread2('Sentinel_Des_LOS.grd');
[x3,y3,z3] = grdread2('Sentinel_Asc_Azi.grd');
[x4,y4,z4] = grdread2('Sentinel_Des_Azi.grd');

Lon = ones(length(y),1)*x;
Lat = flipud(y')*ones(1,length(x));
Inc = [x(2)-x(1),y(2)-y(1)];

%
% Decomposing into E N U 
%

%% Convert2Grd.m
Inv = {Lon,Lat,E,N,U};
Inc = Inc;
Step = 'name'
Convert2Grd(Inv,Inc,Step);
```
##### Go to Terminal

```csh
csh XYZ2GRD.csh
```
You should see *E_name.grd* *N_name.grd* *U_name.grd* in your directory
