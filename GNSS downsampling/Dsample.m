%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Lin,Li-Chieh                                %
%                      Earth and Planetary Sciences                       %
%                  University of California, Riverside                    %
%                              2023.06.17                                 %
%                                                                         %
%                                                                         %
% Down-sampling GNSS or other similar stations due to excessive stations. %
%                                                                         %
% Find stations that are within a defined distance. Remove stations that  %
% have the most stations near them. If they share the same amount of      %
% stations, then take out the one that has the minimum distance among     %
% other stations, meaning that it is closest to its neighbor.             %
%                                                                         %
%                                                                         %
% Input:                                                                  %
% 1. InMat: MxN matrix with local coordinates. First two columns must be  %
% coordinates (X,Y) or (X,Y,Z).                 
% 2. Dist: Stations between stations within this distance will be removed.%
%                                                                         %
% Output:                                                                 %
% 1. OutMat: Remaining stations.                                          %
% 2. OmitInd: The index of the taken-out stations.                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OutMat,NeighborSta] = Dsample(InMat,Dist)
% Call DistMatrix.m to generate distance matrix
DistMat = DistMatrix(InMat);

OutMat = InMat;
NearSta = sum(DistMat<Dist,2)-1;
NeighborSta = length(NearSta(NearSta ~= 0));
MaxStaCount = max(NearSta);
if MaxStaCount ~= 0
    MaxStaInd = find(NearSta == MaxStaCount);
    if length(MaxStaInd) > 1
        Sta = DistMat(MaxStaInd,:);
        ndist = Sta.*(Sta<Dist);
        ndistAvg = sum(ndist,2)/MaxStaCount;
        [~,ndistAvgMinInd] = min(ndistAvg);
        OmitInd = MaxStaInd(ndistAvgMinInd);
        OutMat(OmitInd,:) = [];
    else
        OmitInd = MaxStaInd;
        OutMat(OmitInd,:) = [];
    end
end
end