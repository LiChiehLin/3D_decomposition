%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                         Department of Geography                         %
%                       National Taiwan University                        %
%                              2021.06.01                                 %
%                                                                         %
% Revised: 2021.07.23                                                     %
%          Set additional variable "Datatype" for different data inputs   %
%                                                                         %
% This is to regrid the GRDs from GMTSAR for InSAR3Ddisp.m                %
%                                                                         %
% Note that:                                                              %
% 1. This is an essential step if your GRDs are of different sizes or     %
%    different increments.                                                %
% 2. The data points that fall into the same grid will be averaged        %
%    together.                                                            %
%                                                                         %
% Input:                                                                  %
% 1. InGrd: A cell of Grds (InGrd = {Grd1,Grd2,Grd3...})                  %
% 2. Boundary: A 2-by-2 matrix of the Lower-left and Upper-right point.   %
%              Can take the product from BoundBox.m                       %
% 3. Interval: Increment of the new grid (Interval = [xinc,yinc])         %
%              Can take the product from BoundBox.m                       %
% 4. Num: How many inputs (3 or 4 ....)                                   %
% 5. Datatype: Whether the input is raw or processed(merged)              %
%              0: Raw (Direct product from GMTSAR)                        %
%                 InGrd is a string cell {'grd1','grd2'...}               %
%              1: Processed (Merged Store in .mat file)                   %
%                 InGrd is a big cell {{cell1},{cell2}...}                %
%                                                                         %
%                                                                         %
% Note: When Datatype = 1. Normally is from merging swaths for Sentinel-1 %
%       Otherwise, Datatype is 0 because the sizes and increments from    %
%       Might differ.                                                     %
%                                                                         %
%                                                                         %
% Output:                                                                 %
% 1. Lon: A matrix of the longitude of the new grid                       %
% 2. Lat: A matrix of the latitude of the new grid                        %
% 3. NewGrd: A cell of regridded GRDs (NewGrd = {Grd1,Grd2 ...})          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lon,Lat,NewGrd] = Regrid(InGrd,Boundary,Interval,Num,Datatype)
%% Check inputs
if size(Boundary) ~= [2,2]
    error('Boundary should be a 2-by-2 matrix! ABORT!')
end
if length(Interval) ~= 2
    error('Interval should be a 2-by-1 vector! ABORT!')
end
if length(InGrd) ~= Num
    error('Num should be consistent with the number of input grds! ABORT!')
end

%Declare variables
Bound = Boundary;
Int = Interval;

%Call function MakeGrid.m to generate grid
[grdlon,grdlat,grd]= MakeGrid(Bound,Int);

%Read in all GRDs
outgrd = cell(Num,1);

if Datatype == 0
parfor i = 1:Num
    file = InGrd{i};
    [X,Y,Z] = grdread2(file);
    lonIn = ones(length(Y),1)*X;
    latIn = flipud(Y')*ones(1,length(X));
    ZIn = Z;
    out = ones(length(Y),length(X)); %Store the index of each displacement
    for r = 1:size(latIn,1)
        for c = 1:size(lonIn,2)
            lon = lonIn(r,c);
            lat = latIn(r,c);
            %Conditional statement
            CondLonR = double(lon >= grd(:,2));
            CondLatU = double(lat >= grd(:,3));
            CondLonL = double(lon < grd(:,4));
            CondLatD = double(lat < grd(:,5));
            CondLon = CondLonR == CondLonL;
            CondLat = CondLatU == CondLatD;
            if ~any(CondLon) || ~any(CondLat)
                %If it is out of the grid, then skip
                name = ['Grd:',num2str(i),' ','row:',num2str(r),' ','col:',num2str(c),' ','No match. Out of the grid.'];
                disp(name)
                out(r,c) = nan;
                %continue
            else
                Condmat = [CondLonR,CondLatU,CondLonL,CondLatD];
                Cond = sum(Condmat,2);
                Index = find(Cond == 4);
                name = ['Grd:',num2str(i),' ','row:',num2str(r),' ','col:',num2str(c),' ','In grid! ',num2str(Index)];
                disp(name)
                out(r,c) = Index; %A matrix of new grd index of each data point 
            end
        end
    end
    GrdIndex = ones(size(grd,1),2);
    for j = 1:size(grd,1)
        name = ['Indexing: ',num2str(j)];
        disp(name)
        Index2 = j;
        IndexGrd = find(out == Index2);
        Disp = ZIn(IndexGrd);
        DispMean = mean(Disp(~isnan(Disp)));
        GrdIndex(j,:) = [Index2,DispMean]; 
        %GrdIndex: A matrix of M-by-2. 
        %col1: grd id
        %col2: averaged displacement
    end
    Disp = GrdIndex(:,2);
    %Reshape the grd to make the format of M-by-N
    Regridded = reshape(Disp,size(grdlat,2),size(grdlat,1))';
    outgrd{i} = Regridded;
end

elseif Datatype == 1
parfor i = 1:Num
    file = InGrd{i};
    X = file{1};
    Y = file{2};
    Z = file{3};
    lonIn = X;
    latIn = Y;
    ZIn = Z;
    out = size(ZIn); %Store the index of each displacement
    for r = 1:size(latIn,1)
        for c = 1:size(lonIn,2)
            lon = lonIn(r,c);
            lat = latIn(r,c);
            %Conditional statement
            CondLonR = double(lon >= grd(:,2));
            CondLatU = double(lat >= grd(:,3));
            CondLonL = double(lon < grd(:,4));
            CondLatD = double(lat < grd(:,5));
            CondLon = CondLonR == CondLonL;
            CondLat = CondLatU == CondLatD;
            if ~any(CondLon) || ~any(CondLat)
                %If it is out of the grid, then skip
                name = ['Grd:',num2str(i),' ','row:',num2str(r),' ','col:',num2str(c),' ','No match. Out of the grid.'];
                disp(name)
                out(r,c) = nan;
                %continue
            else
                Condmat = [CondLonR,CondLatU,CondLonL,CondLatD];
                Cond = sum(Condmat,2);
                Index = find(Cond == 4);
                name = ['Grd:',num2str(i),' ','row:',num2str(r),' ','col:',num2str(c),' ','In grid! ',num2str(Index)];
                disp(name)
                out(r,c) = Index; %A matrix of new grd index of each data point 
            end
        end
    end
    GrdIndex = ones(size(grd,1),2);
    for j = 1:size(grd,1)
        name = ['Indexing: ',num2str(j)];
        disp(name)
        Index2 = j;
        IndexGrd = find(out == Index2);
        Disp = ZIn(IndexGrd);
        DispMean = mean(Disp(~isnan(Disp)));
        GrdIndex(j,:) = [Index2,DispMean]; 
        %GrdIndex: A matrix of M-by-2. 
        %col1: grd id
        %col2: averaged displacement
    end
    Disp = GrdIndex(:,2);
    %Reshape the grd to make the format of M-by-N
    Regridded = reshape(Disp,size(grdlat,2),size(grdlat,1))';
    outgrd{i} = Regridded;
end

end
Lon = grdlon;
Lat = grdlat;
NewGrd = outgrd;
end


