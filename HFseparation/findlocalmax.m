%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Lin,Li-Chieh                                %
%                         Department of Geography                         %
%                       National Taiwan University                        %
%                              2021.08.11                                 %
%                                                                         %
%                                                                         %
% Find the local maxima of 1D array                                       %
%                                                                         %
% Input:                                                                  %
% 1. Array: 1D array                                                      %
% 2. Neighbor: Search of adjacent points to find the maxima (number)      %
%              Larger to avoid unwanted maxima resulted from oscillation  %
%              Smaller to include all possible maxima                     %
%                                                                         %
% Output:                                                                 %
% 1. LocalMax: Value of local maxima                                      %
% 2. Index: Index of local maxima                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [LocalMax,Index] = findlocalmax(Array,Neighbor)
Adj = Neighbor;
Arg = Array;

j = 0;
for i = (1+Adj):(length(Array)-Adj)
    Ind = (i-Adj):(i+Adj);
    CenPts = Arg(Ind(Adj+1));
    AdjPts = Arg(setdiff(Ind,ceil(length(Ind)/2)));
    if sum(CenPts > AdjPts) == Adj*2
        j = j + 1;
        Index(j) = Ind(Adj+1);
        LocalMax(j) = CenPts;
    end
end

end