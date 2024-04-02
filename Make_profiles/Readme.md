# Make_profiles

Make profiles based on a *given point*, *length*, *width* and *direction (azimuth)*   

#### Two MatLab codes:
* ProfileLine.m (Generates the profile indices. Sub-routine in `MakeProfile.m` but can be used independently)  
* MakeProfile.m (Samples the input grd with profile indices generated from `ProfileLine.m`)
---
Usage of `MakeProfile.m`  
##### Input variable:  
   * OriginR: Supports one point or a vector of points. This is the Row index.  
   * OriginC: Supports one point or a vector of points. This is the Col index.  
   * Azi: Number. The profile lengthening direction with respect to right.  
     * 0: profile lengthening towards right  
     * 90: profile lengthening towards up
   * Length: Number. The length of the profile. Count by pixels.
   * Width: Number. The width of the profile.
   * Grd: Matrix. The image matrix being sampled.  
##### Output variable:  
   * Rindex: Vector/matrix. Row index/indices  
   * Cindex: Vector/matrix. Col index/indices  
   * Prof: Vector/matrix.  Sampled pixel values  
     * If only input 1 point, then the ouput is a 1xM vector  
     * If input multipple points, then the ouput is a NxM matrix

The input parameters are identical for `ProfileLine.m`. Just omit the input `grd` and the output `Prof`

---
### Example:  
```MatLab
% Read grd file and set up parameters
[x1,y1,z1] = grdread2('Surface_Velocity.grd');

% Starting points
inc = 10;
OriginC = 5:inc:1080;
OriginR = zeros(1,length(OriginC));
OriginR(1) = 380;

slope = tand(35); % Set up starting points along right-up direction
for i = 2:length(OC)
    OriginR(i) = slope*inc+OriginR(i-1);
end
OriginR = round(OriginR);

% Sampling
Azi = 314;
Len = 600;
Width = 0;
[Rindex,Cindex,Prof] = MakeProfile(OriginR,OriginC,Azi,Len,Width,z1);
```
### Visualize:  
```MatLab
figure();hold on
imagesc(z1)
plot(Cindex,Rindex,'r.','MarkerSize',4)
hold off
```
![Example](https://github.com/LiChiehLin/3D_decomposition/blob/8d45c7e17fbb76a73330e99079901795eb7aba88/Figures/Make_Profiles_Example.png)
