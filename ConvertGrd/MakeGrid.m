%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                         Department of Geography                         %
%                       National Taiwan University                        %
%                              2021.06.01                                 %
%                                                                         %
%                                                                         %
% Generate a grid for regridding GMTSAR results for further 3D inversion  %
%                                                                         %
% Input:                                                                  %
% 1. Bound: A 2-by-2 matrix of the Lower-left and Upper-right point.      %
%           Can take the product from BoundBox.m                          %
% 2. Interval: Increment of the new grid (Interval = [xinc,yinc])         %
%           Can take the product from BoundBox.m                          %
%                                                                         %
% Note: This is made to be a sub-routine of Regrid.m                      %
%                                                                         %
% Output:                                                                 %
% 1. grdlon: A matrix of the longitude of the new grid                    %
% 2. grdlat: A matrix of the latitude of the new grid                     %
% 3. grd: A matrix containing the two boundary points                     %
% col1: Grid id                                                           %
% col2: Longitude of the lower-left point                                 %
% col3: Latitude of the lower-left point                                  %
% col4: Longitude of the upper-right point                                %
% col5: Latitude of the upper-right point                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [grdlon,grdlat,grd] = MakeGrid(Bound,Interval)
%%%%Convert interval (km) to decimal degree
%lon 1 degree = 102.4699km
londiff = 102.4699;
%lat 1 degree = 111.3195km  in the Taiwan region
latdiff = 111.3195;
intlon = Interval(1);
intlat = Interval(2);

%%%%%Construct sequence of lon and lat
%%Longitude
lonmin = Bound(1,1);
lonmax = Bound(2,1);
latmin = Bound(1,2);
latmax = Bound(2,2);
lon = lonmin:intlon:lonmax+intlon;
lonout = lon;
lonout(end) = [];
%%Latitude
lat = latmin:intlat:latmax+intlat;
latout = lat;
latout(end) = [];
latout = fliplr(latout);
lat = fliplr(lat);

%%%%Calculate the size of the grid
lenlon = length(lon);
lenlonout = length(lonout);
lenlat = length(lat);
lenlatout = length(latout);
lonmat = ones(lenlatout,1)*lonout;
%latflip = flipud(lat');
latflip = latout';
latmat = latflip*ones(1,lenlonout);
grdlon = lonmat;
grdlat = latmat;


grdsize = (lenlon-1)*(lenlat-1);
id = 1:1:grdsize;
id = id';



%%%%Construct grid matrix
%%%Number each grid cell from bottomleft to upperright
%%col1 : bottomleft longitude
%%col2 : bottomleft latitude
%%col3 : upperright longitude
%%col4 : upperright latitude
col1 = lon(1:lenlon-1)';
col3 = lon(2:lenlon)';
t = ones(lenlon-1,1);
for i = 1:lenlat-1
    col4 = lat(i);
    col2 = lat(i+1);
    tmp1 = [col1.*t col2.*t col3.*t col4.*t];
    grd{i} = tmp1;
end
grd = grd';
grd = cell2mat(grd);
grd = [id grd];

end
