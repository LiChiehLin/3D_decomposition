# InSAR3Ddisp.m (Decompose SAR-derived disp. into 3D disp.)
**Note that:**  
Each input variable should align to the image attributes  

##### Input variables:
- InGrd: 
   * cell variable. Each column in the cell should be a image matrix read by *'grdread2.m'*
- Azimuth: 
   * Satellite inclination degree
   * Ascending and Descending should be the same number (Can also be the azimuth angle of each ground pixel)
   * 12° for Sentinel-1 
   * 8° for ALOS-2
- LookAngle:
   * Satellite looking angle (Can also be the azimuth angle of each ground pixel)
   * 32.9° for Sentinel-1 subswath 1 
   * 38.3° for Sentinel-1 subswath 2 
   * 43.1° for Sentinel-1 subswath 3
- DispType:
   * cell variable. Each input is a string either *'LOS'* or *'Azi'*
- Orbit: 
   * cell variable. Each input is a string either *'Asc'* or *'Des'*
- Num:
   * number. How many images are under converting


##### Output variables:
- Out: cell variable. The order is E N U
- Outcount: matrix variable. Counts how many inputs were used to invert ENU displacement
   * Azimuth displacement is the first digit
   * LOS displacement is the second digit
   * If this pixel is inverted from **2 Azimuth displacement** and **1 LOS displacement**
   * Then this particular Outcount(m,n) will be stored **21**

---
### Example:
Make sure every grd has the same size and increment!  
If not, use `Resamp.csh` to make all grds the same sizes and increments
https://github.com/LiChiehLin/GMTSAR_gadgets/blob/de9d5f5c177badab2a16f0c8c42290506e4789b3/ResampleGrd/Resamp.csh
```MatLab
% Read data
[x,y,AscLOS] = grdread2('Sentinel_Asc_LOS.grd');
[~,~,DesLOS] = grdread2('Sentinel_Des_LOS.grd');
[~,~,AscAzi] = grdread2('Sentinel_Asc_Azi.grd');
[~,~,DesAzi] = grdread2('Sentinel_Des_Azi.grd');


% Enter satellite parameters
AlphaAsc = 12; # or a matrix with the same size as the diplacement matrix, inter-changeable
AlphaDes = 12; # or a matrix with the same size as the diplacement matrix, inter-changeable
ThetaAsc = 43.1; # or a matrix with the same size as the diplacement matrix, inter-changeable
ThetaDes = 32.9; # or a matrix with the same size as the diplacement matrix, inter-changeable

% Prepare inversion
InGrd = {AscLOS,DesLOS,AscAzi,DesAzi};
Azimuth = [AlphaAsc,AlphaDes,AlphaAsc,AlphaDes];
# Or: Azimuth = {AlphaAsc,AlphaDes,AlphaAsc,AlphaDes};
Incidence = [ThetaAsc,ThetaDes,ThetaAsc,ThetaDes];
# Or: Incidence = {ThetaAsc,ThetaDes,ThetaAsc,ThetaDes}
DispType = {'LOS','LOS','Azi','Azi'};
Orbit = {'Asc','Des','Asc','Des'};
Num = 4;

% Decomposing into E N U
[Out,Stderr,count] = InSAR3Ddisp(InGrd,Azimuth,Incidence,DispType,Orbit,Num);

% Extract each component
E = ones(size(AscLOS));
Eerr = ones(size(AscLOS));
N = ones(size(AscLOS));
Nerr = ones(size(AscLOS));
U = ones(size(AscLOS));
Uerr = ones(size(AscLOS));
for i = 1:size(AscLOS,1)
    for j = 1:size(AscLOS,2)
        if length(Out{i,j}) == 1
            E(i,j) = nan;
            Eerr(i,j) = nan;
            N(i,j) = nan;
            Nerr(i,j) = nan;
            U(i,j) = nan;
            Uerr(i,j) = nan;
        else
            E(i,j) = Out{i,j}(1);
            Eerr(i,j) = Stderr{i,j}(1);
            N(i,j) = Out{i,j}(2);
            Nerr(i,j) = Stderr{i,j}(2);
            U(i,j) = Out{i,j}(3);
            Uerr(i,j) = Stderr{i,j}(3);
        end
    end
end

%% Convert to grd
grdwrite2(x,y,E,'E.grd');
grdwrite2(x,y,N,'N.grd');
grdwrite2(x,y,U,'U.grd');
```
