# InSAR3Ddisp.m (Decompose SAR-derived disp. into 3D disp.)
### Note that:  
1. Each input attribute should align with the order of InGrd
2. **Positive value: Away from satellite**
3. **Negative value: Closer to satellite**

---

##### Input variables:
- InGrd: 
   * cell variable. Each column in the cell should be a image matrix read by *'grdread2.m'*
- Weights (Covariance matrix $\Sigma$):
   * Larger the number, lower the weight. Supports three formats:
     * A cell of Grds of same size as InGrd. This is to give different weights on every pixel and every input
     * A Mx1 or 1xM vector assigning the same weight over all pixels for the same input
     * A string "no". Assume uniform weighting for every pixel and every input
- Azimuth: 
   * Satellite inclination degree
   * Ascending and Descending should be the same number (Can also be the azimuth angle for each ground pixel)
     * Ascending Sentinel-1: 348° 
     * Descending Sentinel-1 192° 
- LookAngle:
   * Satellite looking angle (Can also be the incidence angle of each ground pixel)
   * 32.9° for Sentinel-1 subswath 1 
   * 38.3° for Sentinel-1 subswath 2 
   * 43.1° for Sentinel-1 subswath 3
- DispType:
   * cell variable. Each input is a string either *'LOS'* or *'Azi'*
- LookDirection:
   * Looking direction. Looking to the SAR system's left or right


##### Output variables:
- Out: cell variable. The order is E N U
- OutModelVar: cell variable. Contains the model variance of E N U. Model variance is given by: $diag([G^T\Sigma^{-1}G]^{-1})$
- Outcount: matrix variable. Counts how many inputs were used to invert ENU displacement
   * Azimuth displacement is the first digit
   * LOS displacement is the second digit
   * If this pixel is inverted from **2 Azimuth displacement** and **1 LOS displacement**
   * Then this particular Outcount(m,n) will be stored **21**

--- 
### Mathematical expression of this decomposition  
The decomposition is done by inversion using **least squares method** when unknowns < knowns  
The linear relationship connecting LOS and E N U can be formed by:  
<p align="center">
$$d = Gm$$  </p>
<p align="center">
where the $G$ matrix is weighted by $\Sigma^{-1}$. Thus, </p>
<p align="center">
$\Sigma^{-1}d=\Sigma^{-1}Gm$ </p>
<p align="center">
Then, the solution $m$ and model variance $\sigma^{2}$ are given by: </p>
<p align="center">
$m=[G^{T}\Sigma^{-1}G]^{-1}G^{T}\Sigma^{-1}d$ </p>
<p align="center">
$\sigma^{2}=diag([G^{T}\Sigma^{-1}G]^{-1})$ </p>


---
### Example:
Make sure every grd has the same size and increment!  
If not, use `Resamp.csh` to make all grds the same sizes and increments  
Please refer to my another repository: `LiChiehLin/GMTSAR_gadgets/ResampleGrd/Resamp.csh`
```MatLab
% Read data
[x,y,AscLOS] = grdread2('Sentinel_Asc_LOS.grd');
[~,~,DesLOS] = grdread2('Sentinel_Des_LOS.grd');
[~,~,AscAzi] = grdread2('Sentinel_Asc_Azi.grd');
[~,~,DesAzi] = grdread2('Sentinel_Des_Azi.grd');


% Enter satellite parameters
AlphaAsc = 348; # or a matrix with the same size as the diplacement matrix, inter-changeable
AlphaDes = 192; # or a matrix with the same size as the diplacement matrix, inter-changeable
ThetaAsc = 43.1; # or a matrix with the same size as the diplacement matrix, inter-changeable
ThetaDes = 32.9; # or a matrix with the same size as the diplacement matrix, inter-changeable

% Prepare inversion
InGrd = {AscLOS,DesLOS,AscAzi,DesAzi};
Weights = [0.1, 0.1, 1.0, 1.0]; % Larger the number, lower the weight
Azimuth = [AlphaAsc,AlphaDes,AlphaAsc,AlphaDes];
% Or: Azimuth = {AlphaAsc,AlphaDes,AlphaAsc,AlphaDes};
Incidence = [ThetaAsc,ThetaDes,ThetaAsc,ThetaDes];
% Or: Incidence = {ThetaAsc,ThetaDes,ThetaAsc,ThetaDes}
DispType = {'LOS','LOS','Azi','Azi'};
LookDirection = {'r','r','r','r'};

% Decomposing into E N U
[Out,ModelVar,count] = InSAR3Ddisp(InGrd,Weights,Azimuth,Incidence,DispType,LookDirection);

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
            Eerr(i,j) = ModelVar{i,j}(1);
            N(i,j) = Out{i,j}(2);
            Nerr(i,j) = ModelVar{i,j}(2);
            U(i,j) = Out{i,j}(3);
            Uerr(i,j) = ModelVar{i,j}(3);
        end
    end
end

%% Convert to grd
grdwrite2(x,y,E,'E.grd');
grdwrite2(x,y,N,'N.grd');
grdwrite2(x,y,U,'U.grd');
```
