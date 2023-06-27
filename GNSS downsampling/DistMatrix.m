%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Lin,Li-Chieh                                %
%                      Earth and Planetary Sciences                       %
%                  University of California, Riverside                    %
%                              2023.06.17                                 %
%                                                                         %
%                                                                         %
% Calculate distances between each input points. (e.g. GNSS stations)     %
%                                                                         %
% Input:                                                                  %
% 1. InMat: MxN matrix with local coordinates. First two columns must be  %
% coordinates (X,Y) or (X,Y,Z).                                           %                                           
%                                                                         %
% Output:                                                                 %
% 1. OutMat: Distance matrix of each input points                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OutMat = DistMatrix(InMat)
OutMat = zeros(size(InMat,1),size(InMat,1));
for i = 1:size(InMat,1)
    tmpMat = InMat(:,[1,2]);
    tmprow = size(tmpMat,1);
    Now = tmpMat(i,:);
    NowMat = ones(tmprow,1)*Now;
    dist = sqrt(sum((NowMat-tmpMat).^2,2));
    OutMat(i,:) = dist;
end
end