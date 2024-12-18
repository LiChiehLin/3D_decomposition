%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Lin,Li-Chieh                                %
%                     Earth and Planetary Sciences                        %
%                  University of California, Riverside                    %
%                              2024.12.13                                 %
%                                                                         %
% Perform Bit-Plane slicing reconstruction based on the product of        %
% BPslice.m. Given the bit planes, reconstruct the original input gray    %
% image. This aims to get rid of some irrelvant areas, such as noise,     %
% secondary signals in the image                                          %
%                                                                         %
% Input:                                                                  %
% 1. OrigMatrix: A MxN image matrix that was used in BPslice.m. This could%
% also be a different matrix but with the same size. e.g. If you wish to  %
% reconstruct the LOS velocity using only a few bit planes.               %
% 2. BitPlanes: A 8x1 cell Bit-planes from BPslice.m                      %
% 3. BitPlaneInd: An Mx1 array that shows what planes for reconstruction  %
%                                                                         %
% Output:                                                                 %
% 1. Reconstruct: The reconstructed image of the input "OrigMatrix".      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Reconstruct = BPrecon(OrigMatrix,BitPlanes,BitPlaneInd)
disp(strcat('*** Reconstructing based on bit planes',32,num2str(BitPlaneInd)))
N = length(BitPlaneInd);
BPbool = nan(size(OrigMatrix,1)*size(OrigMatrix,2),N);
% Get the indices of each plane
% This might take a little while (a few seconds)
for i = 1:length(BitPlaneInd)
    ind = BitPlaneInd(i);
    tmp = BitPlanes{ind};
    BPbool(:,i) = tmp(:);
end

% Get the combnined indices
BPInd = sum(BPbool(:,1:N),2);
BPInd(BPInd~=0 & ~isnan(BPInd)) = 1;
BPInd(BPInd==0 | isnan(BPInd)) = 0;
BPInd = boolean(BPInd);
Reconstruct = OrigMatrix;
Reconstruct(~BPInd) = nan;

end





