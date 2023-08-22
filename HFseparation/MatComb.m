%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Lin,Li-Chieh                                %
%                         Department of Geography                         %
%                       National Taiwan University                        %
%                              2021.08.17                                 %
%                                                                         %
%                                                                         %
% Combining two matrices with same sizes into one matrix                  %
% A routine for post processing HFsep.m                                   %
%                                                                         %
% Element combining will have three cases:                                %
% 1. Either A or B has value for (i,j), output the value for (i,j)        %
% 2. Neither A nor B has value for (i,j), put NaN                         %
% 3. Both A and B have value for (i,j), take the average of them for (i,j)%
%                                                                         %
% Input:                                                                  %
% 1. MatA: Matrix to be combined                                          %
% 2. MatB: Matrix to be combined                                          %
%                                                                         %
% Output:                                                                 %
% 1. MatComb: Combined matrix                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MatComb = MatComb(MatA,MatB)
if size(MatA) ~= size(MatB)
    error('Matrices of different sizes')
end


[R,C] = size(MatA);
MatComb = zeros(R,C);
for r = 1:R
    for c = 1:C
        PixA = MatA(r,c);
        PixB = MatB(r,c);
        % Case1: A has value B is NaN
        if ~isnan(PixA) && isnan(PixB)
            MatComb(r,c) = PixA;
            
        % Case2: A is NaN B has value
        elseif isnan(PixA) && ~isnan(PixB)
            MatComb(r,c) = PixB;
            
        % Case3: A and B have value
        elseif ~isnan(PixA) && ~isnan(PixB)
            MatComb(r,c) = mean([PixA,PixB]);
            
        % Case4: A and B are NaN
        elseif isnan(PixA) && isnan(PixB)
            MatComb(r,c) = nan;
        end
    end
end






end