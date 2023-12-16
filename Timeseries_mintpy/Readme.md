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