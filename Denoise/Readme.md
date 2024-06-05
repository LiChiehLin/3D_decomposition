# Deno.m (Adaptive Thresholding)

Approximately remove ~5% of data at the current input image matrix.  

1. A self-defined window size slides through the whole image with step=1
2. Calculate the `sum of difference σ` between the center pixel and adjacent pixels of each step
3. Calculate the Cumulative Distribution Function (CDF) of all `σ`
4. Determine `Threshold γ` based on 95% of the CDF
5. Remove pixels that the corresponding `σ` is larger than `γ`

---
##### Input variable:
   * NoiseMat: matrix. Read from *grdread2.m*
   * ws: number. Window size (Has to be an odd number)
##### Output variable:
   * DeNoiseMat: matrix. Denoised image matrix
   * NInd: matrix. The indices of removed pixels
---
#### Please cite the paper if you used the code.  
Lin, L. C. J., Chuang, R. Y., Lu, C. H., Ching, K. E., & Chen, C. L. (2024). *Derivation of 3D Coseismic Displacement Field from Integrated Azimuth and LOS Displacements for the 2018 Hualien Earthquake.* Remote Sensing, 16(7), 1159.

---
### Example:
```MatLab
% Read data and determine window size
[x,y,ToBeRemoved] = grdread2('Displacement.grd');
ws = 3;

% Start program Deno.m
[Removed,NInd] = Deno(ToBeRemoved,ws);
```
![Example](https://github.com/LiChiehLin/3D_decomposition/blob/57d1ca3b4f9f5562579f5f50b04ef2223b09013b/Figures/Denoise_Example.png)
