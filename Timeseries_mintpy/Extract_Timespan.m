%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Lin,Li-Chieh                                %
%                     Earth and Planetary Sciences                        %
%                  University of California, Riverside                    %
%                              2023.12.15                                 %
%                                                                         %
% (2024.01.27) Revision:                                                  %
% Change parameterization in velocity estimation, add checks before start %
%                                                                         %
% Extract surface velocity field in certain time period from timeseries   %
% data made by mintpy (geo_timeseries.h5, geo_timeseries_ramp.h5)         %
%                                                                         %
% Input:                                                                  %
% 1. Inh5file: String. The name of the .h5 data                           %
% 2. StartT: Number. The starting time of your timespan                   %
% 3. EndT: Number. The ending time of your timespan                       %
% Note that the input time format should be 8 digits with the form of     %
% (yyyymmdd)                                                              %
%                                                                         %
%                                                                         %
% Output:                                                                 %
% 1. Lon: Vector. The Longitude (Compatible with grdwrite2)               %
% 2. Lat: Vector. The Latitude (Compatible with grdwrite2)                %
% 3. Out: 3D matrix. The velocity, velocity standard err, incercept err   % 
%                                                                         %
% Usage:                                                                  %
% [Lon,Lat,OutVelo] =                                                     %
% Extract_Timespan('geo_timeseries.h5',20091008,20221117)                 %
%                                                                         %
% Note:                                                                   %
% Use grdwrite2(Lon,Lat,Out(:,:,1),'GRDNAME.grd') to output as .grd file  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lon,Lat,Out] = Extract_Timespan(Inh5file,StartT,EndT)
%% Check input time format
if length(num2str(StartT)) ~= 8 || length(num2str(EndT)) ~= 8
    error('Input time format needs to be 8 digits (yyyymmdd)')
end

%% Read matrix and attributes
file = hdf5info(Inh5file);

% Get dates
epoch = hdf5read(file.GroupHierarchy.Datasets(2));
date = zeros(length(epoch),1);
for i = 1:length(epoch)
    date(i,1) = str2double(epoch(i).Data);
end
disp(strcat('Full length of data:',num2str(date(1)),'~',num2str(date(end))))
disp(strcat('Extracting data between:',num2str(StartT),'~',num2str(EndT)))

% Convert to decimal year
datestr = num2str(date);
yyyy = str2num(datestr(:,1:4));
mm = str2num(datestr(:,5:6));
dd = str2num(datestr(:,7:8));
datedeci = decyear(yyyy,mm,dd);


% Get displacement
Displ = hdf5read(file.GroupHierarchy.Datasets(3));
% Transpose and flipud the matrix to fit with Lon Lat and the format
% grdread2 reads in (flipud when doing imagesc)
Displ = permute(Displ,[2,1,3]); % Transpose
Displ = flipud(Displ); % flipud

% Get Lon Lat
% Note that there's a discrepancy between the coordinate steps and the
% autual ouput. Below fixed this problem
Attr = cell(length(file.GroupHierarchy.Attributes),1);
for i = 1:length(file.GroupHierarchy.Attributes)
    Attr{i,1} = file.GroupHierarchy.Attributes(i).Shortname;
end
Lonstart = str2double(file.GroupHierarchy.Attributes(strcmp(Attr,'X_FIRST')).Value.Data);
Lonstep = str2double(file.GroupHierarchy.Attributes(strcmp(Attr,'X_STEP')).Value.Data);
Lonend = Lonstart+Lonstep*size(Displ,2);
Lonstep = (Lonend - Lonstart)/(size(Displ,2)-1);

Latstart = str2double(file.GroupHierarchy.Attributes(strcmp(Attr,'Y_FIRST')).Value.Data);
Latstep = str2double(file.GroupHierarchy.Attributes(strcmp(Attr,'Y_STEP')).Value.Data);
Latend = Latstart+Latstep*size(Displ,1);
Latstep = (Latend - Latstart)/(size(Displ,1)-1);

Lon = Lonstart:Lonstep:(size(Displ,2)-1)*Lonstep+Lonstart;
Lat = Latstart:Latstep:(size(Displ,1)-1)*Latstep+Latstart;
if (Lat(1) > Lat(end))
    Lat = fliplr(Lat);
end

%% Get displacement for given timespan
Ind = (date >= StartT) & (date <= EndT);
P1 = Displ(:,:,Ind);

%% Estimate velocity using simple linear regression
row = size(P1,1); col = size(P1,2);
P1calc = zeros(size(P1,3),row*col);
for i = 1:size(P1,3)
    tmp = P1(:,:,i);
    P1calc(i,:) = tmp(:);
end
x = datedeci;
G = [x,ones(size(P1,3),1)];
m = G\P1calc;

slope = m(1,:);
intercept = m(2,:);

%% Estimate standard error of slope and intercept
n = size(P1calc,1)*ones(1,size(P1calc,2));
Sx = sum(x)*ones(1,size(P1calc,2)); Sy = sum(P1calc,1);
Sxx = sum(x.^2)*ones(1,size(P1calc,2)); Syy = sum(P1calc.^2,1);
Se2 = (1./(n.*(n-2)).*(n.*Syy-Sy.^2-slope.^2.*(n.*Sxx-Sx.^2)));
Sb2 = (n.*Se2)./(n.*Sxx - Sx.^2);
Sa2 = Sb2.*(1./n).*Sxx;
aerr = Tlookup(0.95,n(1)-2).*sqrt(Sa2);
berr = Tlookup(0.95,n(1)-2).*sqrt(Sb2);


OutVelo = reshape(slope,row,col);
OutSlopeErr = reshape(berr,row,col);
OutInterceptErr = reshape(aerr,row,col);


Out(:,:,1) = OutVelo;
Out(:,:,2) = OutSlopeErr;
Out(:,:,3) = OutInterceptErr;
end




















