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
% (2024.03.14) Update:                                                    %
% 1. Change from hdf5read to h5read to make it more efficient             %
% 2. Change how it reads data matrices to not take up too much memory     %
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

%% Read dates
file = h5info(Inh5file);

epoch = h5read(Inh5file,'/date');
for i = 1:length(epoch)
    date(i,1) = str2double(epoch(i));
end
disp(date)

% Convert to decimal year
datestr = num2str(date);
yyyy = str2num(datestr(:,1:4));
mm = str2num(datestr(:,5:6));
dd = str2num(datestr(:,7:8));
datedeci = decyear(yyyy,mm,dd);

%% Get displacement for given timespan
Ind = (date >= StartT) & (date <= EndT);
Sind = find(Ind); Sind = Sind(1);
Eind = find(Ind); Eind = Eind(end);
Count = length(find(Ind));
disp(strcat('Extracting data between:',32,num2str(date(Sind)),'~',num2str(date(Eind))))

% Get displacement
dim = file.Datasets(3).Dataspace.Size;
Displ = h5read(Inh5file,'/timeseries',[1 1 Sind],[dim(1) dim(2) Count]);

% Transpose and flipud the matrix to fit with Lon Lat and the format
% grdread2 reads in (flipud when doing imagesc)
Displ = permute(Displ,[2,1,3]); % Transpose
Displ = flipud(Displ); % flipud

% Get Lon Lat
% Note that there's a discrepancy between the coordinate steps and the
% autual ouput. Below fixed this problem
Attr = cell(length(file.Attributes),1);
for i = 1:length(file.Attributes)
    Attr{i,1} = file.Attributes(i).Name;
end
Lonstart = str2double(file.Attributes(strcmp(Attr,'X_FIRST')).Value);
Lonstep = str2double(file.Attributes(strcmp(Attr,'X_STEP')).Value);
Lonend = Lonstart+Lonstep*size(Displ,2);
Lonstep = (Lonend - Lonstart)/(size(Displ,2)-1);

Latstart = str2double(file.Attributes(strcmp(Attr,'Y_FIRST')).Value);
Latstep = str2double(file.Attributes(strcmp(Attr,'Y_STEP')).Value);
Latend = Latstart+Latstep*size(Displ,1);
Latstep = (Latend - Latstart)/(size(Displ,1)-1);

Lon = Lonstart:Lonstep:(size(Displ,2)-1)*Lonstep+Lonstart;
Lat = Latstart:Latstep:(size(Displ,1)-1)*Latstep+Latstart;
if (Lat(1) > Lat(end))
    Lat = fliplr(Lat);
end


%% Estimate velocity using simple linear regression
x = datedeci(Ind);
row = size(Displ,1); col = size(Displ,2);
P1calc = zeros(size(Displ,3),row*col);
for i = 1:size(Displ,3)
    tmp = Displ(:,:,i);
    P1calc(i,:) = tmp(:);
end
G = [x,ones(size(Displ,3),1)];
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




















