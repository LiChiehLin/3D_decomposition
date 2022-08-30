%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                         Department of Geography                         %
%                       National Taiwan University                        %
%                              2021.08.09                                 %
%                                                                         %
% Ver.2: 2021.08.18                                                       %
% Add noise index output (NInd)                                           %
%                                                                         %
%                                                                         %
% (Adaptive Thresholding)                                                 %
% Remove noises from original matrix                                      %
% Set unwanted pixels NaN if the difference between the center pixel and  %
% the adjacent pixels are too big with respect to a relative noise-free   %
% window (i.e. 3x3)                                                       %
%                                                                         %
% Input:                                                                  %
% 1. NoiseMat: Matrix with pixels that contain noises                     %
% 2. ws: Window under calculation (has to be an odd noumber)              %
%                                                                         %
% Output:                                                                 %
% 1. DeNoiseMat: Matrix that has been de-noised                           %
% 2. NInd: NaN indices                                                    %
%                                                                         %
% Note that: This code approximately set 5% of NaN among the whole data   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DeNoiseMat,NInd] = Deno(NoiseMat,ws)
[R,C] = size(NoiseMat);
Grd = NoiseMat;

% Determine the edge of the matrix with respect to the window size
Edge = floor(ws/2);
% Construct boundary row and column
for i = 1:Edge
    RowFirst(i) = i;
    RowEnd(i) = R - (i-1);
    ColFirst(i) = i;
    ColEnd(i) = C - (i-1);
end
EdgeRow = [RowFirst,fliplr(RowEnd)];
EdgeCol = [ColFirst,fliplr(ColEnd)];

% Determine the relative noise-free window
DiffAvg = zeros(length((1+Edge):(R-Edge))*length((1+Edge):(C-Edge)),1);
i = 0;
for r = (1+Edge):(R-Edge)
    for c = (1+Edge):(C-Edge)
        i = i + 1;
        IndexR = r-Edge:r+Edge;
        IndexC = c-Edge:c+Edge;
        
        CenterPix = Grd(r,c);
        MatCal = Grd(IndexR,IndexC);
        MatCal(isnan(MatCal)) = [];
        Diff = abs(CenterPix - MatCal);
        DiffAvg(i) = mean(Diff(:));
    end
end
DiffAvg(isnan(DiffAvg)) = [];

% Compute CDF of DiffAvg to determine the Tolerance

  % Slow when using the following three
%h = histogram(DiffAvg);
%Counts = h.BinCounts;
%Bins = h.BinEdges;
[Counts,Bins] = histcounts(DiffAvg); % Use this otherwise
Bins(1) = [];

CDF = cumsum(Counts)/sum(Counts);

% Tolerance is the difference of 95%
Index = abs(CDF - 0.95);
[~,Index] = min(Index);
Tol = Bins(Index);

out = zeros(size(Grd));
i = 0;
% Deal with inner pixels first
for r = (1 + Edge):(R - Edge)
    for c = (1 + Edge):(C - Edge)
        % Window that the Pix is centered on
        IndexR = r-Edge:r+Edge;
        IndexC = c-Edge:c+Edge;
        CenterPix = Grd(r,c);
        MatCal = Grd(IndexR,IndexC);
        MatCal(isnan(MatCal)) = [];
        Diff = abs(CenterPix - MatCal);
        DiffAvg = mean(Diff(:));
        if DiffAvg > Tol
            i = i + 1;
            % Difference larger than tolerance
            CenterPix = nan;
            out(r,c) = CenterPix;
            outind(i) = sub2ind(size(out),r,c);
        else
            % Difference smaller than tolerance
            out(r,c) = CenterPix;
        end
    end 
end

% Deal with boundary pixels
for r = 1:R
    if any(r == EdgeRow)
        for c = 1:C
            if find(r == EdgeRow) <= (length(EdgeRow)/2)
                IndexR = r-Edge:r+Edge;
                IndexR(IndexR <= 0) = [];
                if find(c == EdgeCol) <= (length(EdgeCol)/2)
                    % Upperleft corner
                    IndexC = c-Edge:c+Edge;
                    IndexC(IndexC <= 0) = [];
                    CenterPix = Grd(r,c);
                    MatCal = Grd(IndexR,IndexC);
                    MatCal(isnan(MatCal)) = [];
                    Diff = abs(CenterPix - MatCal);
                    DiffAvg = mean(Diff(:));
                    if DiffAvg > Tol
                        i = i + 1;
                        % Difference larger than tolerance
                        CenterPix = nan;
                        out(r,c) = CenterPix;
                        outind(i) = sub2ind(size(out),r,c);
                    else
                        % Difference smaller than tolerance
                        out(r,c) = CenterPix;
                    end
                    
                elseif find(c == EdgeCol) > (length(EdgeCol)/2)
                    % Upperright corner
                    IndexC = c-Edge:c+Edge;
                    IndexC(IndexC > C) = [];
                    CenterPix = Grd(r,c);
                    MatCal = Grd(IndexR,IndexC);
                    MatCal(isnan(MatCal)) = [];
                    Diff = abs(CenterPix - MatCal);
                    DiffAvg = mean(Diff(:));
                    if DiffAvg > Tol
                        i = i + 1;
                        % Difference larger than tolerance
                        CenterPix = nan;
                        out(r,c) = CenterPix;
                        outind(i) = sub2ind(size(out),r,c);
                    else
                        % Difference smaller than tolerance
                        out(r,c) = CenterPix;
                    end
                    
                else
                    % Upper boundary
                    IndexC = c-Edge:c+Edge;
                    CenterPix = Grd(r,c);
                    MatCal = Grd(IndexR,IndexC);
                    MatCal(isnan(MatCal)) = [];
                    Diff = abs(CenterPix - MatCal);
                    DiffAvg = mean(Diff(:));
                    if DiffAvg > Tol
                        i = i + 1;
                        % Difference larger than tolerance
                        CenterPix = nan;
                        out(r,c) = CenterPix;
                        outind(i) = sub2ind(size(out),r,c);
                    else
                        % Difference smaller than tolerance
                        out(r,c) = CenterPix;
                    end
                    
                end
            else
                IndexR = r-Edge:r+Edge;
                IndexR(IndexR > R) = [];
                if find(c == EdgeCol) <= (length(EdgeCol)/2)
                    % Lowerleft corner
                    IndexC = c-Edge:c+Edge;
                    IndexC(IndexC <= 0) = [];
                    CenterPix = Grd(r,c);
                    MatCal = Grd(IndexR,IndexC);
                    MatCal(isnan(MatCal)) = [];
                    Diff = abs(CenterPix - MatCal);
                    DiffAvg = mean(Diff(:));
                    if DiffAvg > Tol
                        i = i + 1;
                        % Difference larger than tolerance
                        CenterPix = nan;
                        out(r,c) = CenterPix;
                        outind(i) = sub2ind(size(out),r,c);
                    else
                        % Difference smaller than tolerance
                        out(r,c) = CenterPix;
                    end
                    
                elseif find(c == EdgeCol) > (length(EdgeCol)/2)
                    % Lowerright corner
                    IndexC = c-Edge:c+Edge;
                    IndexC(IndexC > C) = [];
                    CenterPix = Grd(r,c);
                    MatCal = Grd(IndexR,IndexC);
                    MatCal(isnan(MatCal)) = [];
                    Diff = abs(CenterPix - MatCal);
                    DiffAvg = mean(Diff(:));
                    if DiffAvg > Tol
                        i = i + 1;
                        % Difference larger than tolerance
                        CenterPix = nan;
                        out(r,c) = CenterPix;
                        outind(i) = sub2ind(size(out),r,c);
                    else
                        % Difference smaller than tolerance
                        out(r,c) = CenterPix;
                    end
                    
                else
                    % Lower boundary
                    IndexC = c-Edge:c+Edge;
                    CenterPix = Grd(r,c);
                    MatCal = Grd(IndexR,IndexC);
                    MatCal(isnan(MatCal)) = [];
                    Diff = abs(CenterPix - MatCal);
                    DiffAvg = mean(Diff(:));
                    if DiffAvg > Tol
                        i = i + 1;
                        % Difference larger than tolerance
                        CenterPix = nan;
                        out(r,c) = CenterPix;
                        outind(i) = sub2ind(size(out),r,c);
                    else
                        % Difference smaller than tolerance
                        out(r,c) = CenterPix;
                    end
                    
                end
            end
            
            
        end
    else
        
        for c = EdgeCol
            IndexR = r-Edge:r+Edge;
            if find(c == EdgeCol) <= (length(EdgeCol)/2)
                % Left
                IndexC = c-Edge:c+Edge;
                IndexC(IndexC <= 0) = [];
                CenterPix = Grd(r,c);
                MatCal = Grd(IndexR,IndexC);
                MatCal(isnan(MatCal)) = [];
                Diff = abs(CenterPix - MatCal);
                DiffAvg = mean(Diff(:));
                if DiffAvg > Tol
                    i = i + 1;
                    % Difference larger than tolerance
                    CenterPix = nan;
                    out(r,c) = CenterPix;
                    outind(i) = sub2ind(size(out),r,c);
                else
                    % Difference smaller than tolerance
                    out(r,c) = CenterPix;
                end
            else
                % Right
                IndexC = c-Edge:c+Edge;
                IndexC(IndexC > C) = [];
                CenterPix = Grd(r,c);
                MatCal = Grd(IndexR,IndexC);
                MatCal(isnan(MatCal)) = [];
                Diff = abs(CenterPix - MatCal);
                DiffAvg = mean(Diff(:));
                if DiffAvg > Tol
                    i = i + 1;
                    % Difference larger than tolerance
                    CenterPix = nan;
                    out(r,c) = CenterPix;
                    outind(i) = sub2ind(size(out),r,c);
                else
                    % Difference smaller than tolerance
                    out(r,c) = CenterPix;
                end
            end
        end
    end
end
DeNoiseMat = out;
NInd = outind;
end



