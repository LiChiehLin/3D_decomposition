# 3D_decomposition
This repository contains a few MatLab codes for converting SAR-derived displacement into true 3D displacement (EW, NS and UD directions) and some useful tools to handle displacement data  

These codes were developed under MatLab 2015b  
  
Put every code in one directory and add path to your *main* script.  
```MatLab
addpath(genpath('DIRECTORY OF MATLAB CODES/'))
```
 
---
Use grdread2.m to read netCDF grds
> The code for reading netCDF grds is from Kelsey Jordahl (2010). grdread2 (https://www.mathworks.com/matlabcentral/fileexchange/25683-grdread2), MATLAB Central File Exchange. Retrieved August 30, 2022.  

---
Use grdwrite2.m to write netCDF grds
> The code for writing netCDF grds is from Kelsey Jordahl (2011). grdwrite2 (https://www.mathworks.com/matlabcentral/fileexchange/26290-grdwrite2), MATLAB Central File Exchange. Retrieved October 20, 2023.

### Below lists the codes in each directories
- 3D decomposition: 
   * InSAR3Ddisp.m  
   * Convert2Grd.m
- Denoise:  
   * Deno.m **(py)**  
- Sptial fitering:  
   * filtsp.m  
- Near neighbor:
  * nearneighbor.m **(py)**  
- GNSS downsampling:
   * DistMatrix.m **(py)**
   * Dsample.m **(py)**
- HFseparation:
  * HFsep.m
  * MatComb.m
  * findlocalmax.m
  * findlocalmin.m
- Timeseries_mintpy:
  * ~~Extract_Timespan.m~~ `(Moved to LiChiehLin/MintPy_gadgets/Extract_velocity)`
  * ~~Tlookup.m~~ `(Moved to LiChiehLin/MintPy_gadgets/Extract_velocity)`
- Expand_extent
  * Expand_extent.m
- Make_profiles
  * ProfileLine.m
  * MakeProfile.m
- Bit-plane slicing
  * BPslice.m
  * BPrecon.m
  * BPshow.m

---
## Python
Python version of the above codes are updating!  
Those have been included in `DispLCL.py` will be marked **(py)** in the above code list  
Include `DispLCL.py` file in your directory, import the function and execute the function 

### Use `DispLCL.py` as follows:
```python
## Import function
import DispLCL
## Execute the function
Output = DispLCL.DistMatrix(input_variables)
```

### For handling netCDF grd file  
```python
import netCDF4 as nc
import matplotlib.pyplot as plt
# Read in and convert to numpy array
f = nc.Dataset('Your_grd')
x = np.asarray(f.variables['x'])
y = np.asarray(f.variables['y'])
z = np.asarray(f.variables['z'])

# Show figure
plt.matshow(z)
```
### For handling .txt file  
```python
data = np.loadtxt('Your_txt',dtype='f')
```

### Required libraries
- numpy
- math
- warnings
