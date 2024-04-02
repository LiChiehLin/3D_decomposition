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
%                                                                         %
% Output:                                                                 %
% 1. Rindex: Row index/indices                                            %
%    If one point, then it's 1xM vector                                   %
%    If multiple points, then it's NxM matrix                             %
% 2. Cindex: Col index/indices                                            %
%    Same as Rindex                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Rindex,Cindex] = ProfileLine(OriginR,OriginC,Azi,Length,Width)
% Make some checks first
if length(OriginR) ~= length(OriginC)
    error('OriginX and OriginY should have the same length')
end

N = length(OriginR);
disp(strcat('Generating',32,num2str(N),32,'profiles...'))

Rindex = zeros(N,Length*(Width*2+1));
Cindex = zeros(N,Length*(Width*2+1));
for i = 1:N
    Rinit = round(OriginR(i)-Width*abs(cosd(Azi)));
    Rend = round(OriginR(i)+Width*abs(cosd(Azi)));
    Cinit = round(OriginC(i)-Width*abs(sind(Azi)));
    Cend = round(OriginC(i)+Width*abs(sind(Azi)));
    R = Rinit:1:Rend;
    C = Cinit:1:Cend;
    if length(R) > length(C)
        R = R';
        C = [C;C(end)*ones(length(R)-length(C),1)];
    elseif length(R) < length(C)
        R = [R,R(end)*ones(length(C)-length(R),1)];
        C = C';
    end

    RowEnd = round(R + Length*sind(Azi));
    ColEnd = round(C + Length*cosd(Azi));
    Rseqtmp = zeros(length(R),Length);
    Cseqtmp = zeros(length(C),Length);
    for j = 1:length(R)
        if R(j) == RowEnd(j)
            Rseq = R(j)*ones(1,Length);
            Cseq = C(j):(ColEnd(j)-C(j))/(Length-1):ColEnd(j);
        elseif C(j) == ColEnd(j)
            Rseq = R(j):(RowEnd(j)-R(j))/(Length-1):RowEnd(j);
            Cseq = C(j)*ones(1,Length);
        else
            Rseq = R(j):(RowEnd(j)-R(j))/(Length-1):RowEnd(j);
            Cseq = C(j):(ColEnd(j)-C(j))/(Length-1):ColEnd(j);
        end
        Rseqtmp(j,:) = round(Rseq);
        Cseqtmp(j,:) = round(Cseq);

    end

    Rindex(i,:) = reshape(Rseqtmp',1,[]);
    Cindex(i,:) = reshape(Cseqtmp',1,[]);

end

end



