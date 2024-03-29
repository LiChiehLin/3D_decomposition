%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Lin,Li-Chieh                                %
%                      Earth and Planetary Sciences                       %
%                  University of California, Riverside                    %
%                              2024.01.09                                 %
%                                                                         %
% Updates:                                                                %
% (2024.02.16) Accomodate grds with slight difference x and y increments  %
% (2024.02.19) Display some stats and progress                            %
%                                                                         %
%                                                                         %
% Expand the input grds to the same spatial extent. This is useful when   %
% doing image comparisons like differencing or adding.                    %
%                                                                         %
% e.g.                                                                    %
% grd1 has boudaries [-122 -121 30 31]                                    %
% grd2 has boudaries [-122.5 -121.5 30.5 31.5]                            %
% The output boudaries will be [-122.5 -121 30 31.5]                      %
%                                                                         %
% Note that: the input grds need to have the same X and Y increments,     %
% otherwise it would cause mis-registration                               %
%-------------------------------------------------------------------------%
% Input:                                                                  %
% 1. InMat: cell. each row containing the x,y,z of the input grds         %
%    col1: x (Lon)                                                        %
%    col2: y (Lat)                                                        %
%    col3: z (displacement/velocity)                                      %
% Note that the format has to be the format read by grdread2.m            %
%                                                                         %
% Output:                                                                 %
% 1. OutMat: cell containing the matrices that have the spatial extent    %
%    row1: x (Expanded Lon)                                               %
%    row2: y (Expanded Lat)                                               %
%    row3~end: Expanded grds                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OutMat = Expand_extent(InMat)
if ~iscell(InMat)
    error('Input should be a cell containing the grds')
end

% Get four boundaries and the mean x and y increments
xmintmp = zeros(size(InMat,1),1);
xmaxtmp = zeros(size(InMat,1),1);
ymintmp = zeros(size(InMat,1),1);
ymaxtmp = zeros(size(InMat,1),1);
xinctmp = zeros(size(InMat,1),1);
yinctmp = zeros(size(InMat,1),1);
for i = 1:size(InMat,1)
    xtmp = InMat{i,1}; ytmp = InMat{i,2};
    xmintmp(i) = min(xtmp); xmaxtmp(i) = max(xtmp);
    ymintmp(i) = min(ytmp); ymaxtmp(i) = max(ytmp);
    xinctmp(i) = abs(xtmp(2) - xtmp(1));
    yinctmp(i) = abs(ytmp(2) - ytmp(1));
    disp(strcat('Original stats for matrix:',num2str(i)))
    disp(strcat('W E S N: ',num2str(min(xtmp)),32,num2str(max(xtmp)),32,num2str(min(ytmp)),32,num2str(max(ytmp))))
    disp(strcat('xinc & yinc: ',num2str(abs(xtmp(2) - xtmp(1))),32,num2str(abs(ytmp(2) - ytmp(1)))))
    disp('--------------------------------------------')
end
xmin = min(xmintmp); xmax = max(xmaxtmp);
ymin = min(ymintmp); ymax = max(ymaxtmp);
xinc = mean(xinctmp); yinc = mean(yinctmp);
disp('New bounding box and increments:')
disp(strcat('W E S N:',num2str(xmin),32,num2str(xmax),32,num2str(ymin),32,num2str(ymax)))
disp(strcat('xinc & yinc:',num2str(xinc),32,num2str(yinc)))
disp('--------------------------------------------')

% Construct Lon and Lat matrices
Lon = xmin:xinc:xmax;
Lat = ymin:yinc:ymax;
LonMat = ones(length(Lat),1)*Lon;
LatMat = Lat'*ones(1,length(Lon));

OutMat{1,1} = Lon;
OutMat{2,1} = Lat;
% Expand original x y to make indices for filling in the grids
for i = 1:size(InMat,1)
    disp(strcat('Expanding matrix:',num2str(i)))
    z = InMat{i,3};
    % Zero-padding for original x and y
    x = InMat{i,1}; y = InMat{i,2};
    [~,xalignind] = min(abs(Lon - x(1)));
    % Force x has the same length as Lon
    if xalignind > (length(Lon) - length(x) + 1)
        xalignind = length(Lon) - length(x) + 1;
    end
    xalign = [nan(1,xalignind-1),x,nan(1,length(Lon)-length(x)-(xalignind-1))];

    [~,yalignind] = min(abs(Lat - y(1)));
    % Force y has the same length as Lat
    if yalignind > (length(Lat) - length(y) + 1)
        yalignind = length(Lat) - length(y) + 1;
    end
    yalign = [nan(1,yalignind-1),y,nan(1,length(Lat)-length(y)-(yalignind-1))];
    
    % Construct the zero-padded to matrices
    xMat = ones(length(yalign),1)*xalign;
    yMat = yalign'*ones(1,length(xalign));

    % Make indices to fill in the grids
    Londiff = LonMat - xMat;
    Latdiff = LatMat - yMat;
    Lonind = ~isnan(Londiff);
    Latind = ~isnan(Latdiff);
    fillind = Lonind & Latind;

    zexpand = nan(size(LonMat));
    zexpand(fillind) = z;
    OutMat{i+2,1} = zexpand;
end

end


