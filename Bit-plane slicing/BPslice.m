%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Lin,Li-Chieh                                %
%                     Earth and Planetary Sciences                        %
%                  University of California, Riverside                    %
%                              2024.12.13                                 %
%                                                                         %
% Perform Bit-Plane slicing on a 0~255 gray scale image                   %
% I followed the way this website does:                                   %
% https://www.geeksforgeeks.org/extract-bit-planes-image-matlab/          %
%                                                                         %
% Input:                                                                  %
% 1. GrayMatrix: A MxN image matrix that the value is in 0~255. If not see%
%    Matlab function "rescale" to see how to convert to [0,1] and further %
%    multiply this by 255.                                                %
%                                                                         %
% Output:                                                                 %
% 1. BitPlanes: A 8x1 cell contains the Bit planes. 1 is the LSB, 8 is the%
%    MSB                                                                  %
% 2. Reconstruct: The reconstructed image based on the bit planes         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BitPlanes,Reconstruct] = BPslice(GrayMatrix)
BitPlanes = cell(8,1);

% Bit-plane slicing
b1 = mod(GrayMatrix,2);
b2 = mod(floor(GrayMatrix/2),2);
b3 = mod(floor(GrayMatrix/4),2);
b4 = mod(floor(GrayMatrix/8),2);
b5 = mod(floor(GrayMatrix/16),2);
b6 = mod(floor(GrayMatrix/32),2);
b7 = mod(floor(GrayMatrix/64),2);
b8 = mod(floor(GrayMatrix/128),2);


BitPlanes{1} = b1;
BitPlanes{2} = b2;
BitPlanes{3} = b3;
BitPlanes{4} = b4;
BitPlanes{5} = b5;
BitPlanes{6} = b6;
BitPlanes{7} = b7;
BitPlanes{8} = b8;

% Reconstructing the image
Reconstruct = (2 * (2 * (2 * (2 * (2 * (2 * (2 * b8 + b7) + b6) + b5) + b4) + b3) + b2) + b1); 
end