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
### Example:
```MatLab
% Read data and determine window size
[x,y,ToBeRemoved] = grdread2('Displacement.grd');
ws = 3;

% Start program Deno.m
[Removed,NInd] = Deno(ToBeRemoved,ws);
```
