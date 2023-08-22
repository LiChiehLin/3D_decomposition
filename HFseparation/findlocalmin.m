%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Lin,Li-Chieh                                %
%                         Department of Geography                         %
%                       National Taiwan University                        %
%                              2021.08.11                                 %
%                                                                         %
%                                                                         %
% Find the local minima of 1D array                                       %
%                                                                         %
% Input:                                                                  %
% 1. Array: 1D array                                                      %
% 2. Neighbor: Search of adjacent points to find the minima (number)      %
%              Larger to avoid unwanted minima resulted from oscillation  %
%              Smaller to include all possible minima                     %
%                                                                         %
% Output:                                                                 %
% 1. LocalMin: Value of local minima                                      %
% 2. Index: Index of local minima                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [LocalMin,Index] = findlocalmin(Array,Neighbor)
Adj = Neighbor;
Arg = Array;

j = 0;
for i = (1+Adj):(length(Array)-Adj)
    Ind = (i-Adj):(i+Adj);
    CenPts = Arg(Ind(Adj+1));
    AdjPts = Arg(setdiff(Ind,ceil(length(Ind)/2)));
    if sum(CenPts < AdjPts) == Adj*2
        j = j + 1;
        Index(j) = Ind(Adj+1);
        LocalMin(j) = CenPts;
    end
end

end