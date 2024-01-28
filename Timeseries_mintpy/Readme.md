# Extract_timespan.m 

To create surface velocity of certain studied period  
From timeseries(displacement) obtained from `mintpy`

##### `Tlookup.m` is a sub-routine in `Extract_timespan.m`. Must be included in path.

---
##### Input variable:
   * Inh5file: string. The name of the .h5 file
   * StartT: number. The starting time of the timespan
   * EndT: number. The ending time of the timespan
##### Output variable:
   * Lon: Vector. Longitude (Can be used in `grdwrite2.m`)
   * Lat: Vector. Latitude (can be used in `grdwrite2.m`)
   * Out: 3D matrix. Each page is as follows:  
     * Velocity
     * Velocity standard error
     * Intercept standar error
---
### Example:
```MatLab
h5 = 'geo_timeseries.h5'
StartT = 20151008;
EndT = 20221117;

% Start program Extract_Timespan.m
[Lon,Lat,Out] = Extract_Timespan(h5,StartT,EndT);
```
#### Output as grd
```MatLab
grdwrite2(Lon,Lat,Out(:,:,1))
```
---
Note that since `MintPy` is built on Python, there exists a little difference due to how Matlab and Python treat decimal values differently.  
The residual is the difference of the above velocities. The unit is mm/yr, so the difference should be negligible.

![Example](https://github.com/LiChiehLin/3D_decomposition/blob/3d79d897e99a70c1cc778293bf0c9dbc7c2f382a/Figures/Extract_Timespan_example.png)
