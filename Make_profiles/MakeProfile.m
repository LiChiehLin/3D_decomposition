%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Lin,Li-Chieh                                %
%                      Earth and Planetary Sciences                       %
%                  University of California, Riverside                    %
%                              2023.03.21                                 %
%                                                                         %
% Make profile indices for image matrices                                 %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. OriginR: Support either one point or multiple points. Row index      %
%    1 point: OriginR = 200                                               %
%    2 or more points: OriginR = [200,300,400,500]                        %
% 2. OriginC: Support either one point or multiple points. Col index      %
%    1 point: OriginR = 300                                               %
%    2 or more points: OriginR = [400,500,600,700]                        %
% Note that OriginR and OriginC should have the same length               %
% 3. Azi: Profile direction. Angle w.r.t. right and counter-clockwise     %
%    Right: 0 degree                                                      %
%    Up: 90 degree                                                        %
% 4. Length: The length of the profile (Count by pixels)                  %
% 5. Width: The width of the profile. 0 is one line, 1 is three lines etc.%
% 6. Grd: The image matrix being sampled                                  %
%                                                                         %
% Output:                                                                 %
% 1. Rindex: Row index/indices                                            %
%    If one point, then it's 1xM vector                                   %
%    If multiple points, then it's NxM matrix                             %
% 2. Cindex: Col index/indices                                            %
%    Same as Rindex                                                       %
% 3. Prof: The pixel values along the profile                             %
%    Same as Rindex                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rindex,Cindex,Prof] = MakeProfile(OriginR,OriginC,Azi,Length,Width,Grd)
% Make profile lines
[Rindex, Cindex] = ProfileLine(OriginR,OriginC,Azi,Length,Width);

% Check if any profile line is longer than the image matrix
Row = size(Grd,1);
Col = size(Grd,2);
Rindex(Rindex <= 0) = nan;
Rindex(Rindex > Row) = nan;
Cindex(Cindex <= 0) = nan;
Cindex(Cindex > Col) = nan;

% Convert to linear indexing
N = size(Rindex,1);
M = size(Rindex,2);
Prof = zeros(N,M);
for n = 1:N
    R = Rindex(n,:);
    C = Cindex(n,:);
    % Ensure no nan indices and store nan if index is nan
    lind = sub2ind(size(Grd),R,C);
    for i = 1:length(lind)
        if isnan(lind(i))
            Prof(n,i) = nan;
        else
            Prof(n,i) = Grd(lind(i));
        end
    end

end

end






