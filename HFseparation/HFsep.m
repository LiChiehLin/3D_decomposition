%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Lin,Li-Chieh                                %
%                         Department of Geography                         %
%                       National Taiwan University                        %
%                              2021.08.11                                 %
%                                                                         %
%                                                                         %
% Separates hanging wall and foot wall into two matrices                  %
% Finds the peak value of a subset of the input matrix by grid searching  %
% given the strike of the fault                                           %
%                                                                         %
% This function takes the characteristic of the different deformation     %
% pattern from hanging wall and foot wall to do the separation. So        %
% preliminary de-noising might be needed before-hand                      %
%                                                                         %
% Input:                                                                  %
% 1. InGrd: Matrix that needs to be separated                             %
% 2. Strike: Strike direction of the fault. Or the edge direction of the  %
%            input matrix. Degree with respect to north.                  %
%            i.e. North-South fault: 0 degrees                            %
%                 Northeast fault: 45 degrees                             %
%                 East-West fault: 90 degrees                             %
% 3. PtoN: Matrix value going from positive to negative or not            %
%          1: Yes, 0: No                                                  %
%          Function searches perpendicular to fault strike direction, so  %
%          Strike = 0: Searches from left to right                        %
%                      Put 1 if left is positive, 0 if left is negative   %
%          Strike = 90: Searches from up to down                          %
%                      Put 1 if up is positive, 0 if up is negative       %
% 4. Overlap: How much overlapping is there between the two walls         %
%             zero overlapping is not recommended                         %
%                                                                         %
%                                                                         %
% Output:                                                                 %
% 1. SideA: one wall of the matrix                                        %
% 2. SideB: the other wall of the matrix                                  %
% 3. SideAIndex: matrix index of SideA                                    %
% 4. SideBIndex: matrix index of SideB                                    %
%                                                                         %
% Note that: 1. Size of output matrices are identical with input for      %
%               convenient use                                            %
%            2. Output indices are linear indexing                        %
%               SEE ALSO ind2sub/sub2ind                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SideA,SideB,SideAIndex,SideBIndex] = HFsep(InGrd,Strike,PtoN,Overlap)
if Strike >= 0 && Strike <= 90
    Str = Strike;
    Grd = InGrd;
    fl = 0;
elseif Strike > 90 && Strike <= 180
    Str = 180 - Strike;
    Grd = fliplr(InGrd);
    fl = 1;
elseif Strike > 180 && Strike <= 270 
    Str = Strike - 180;
    Grd = InGrd;
    fl = 0;
elseif Strike > 270 && Strike <= 360
    Str = 360 - Strike;
    Grd = fliplr(InGrd);
    fl = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Some notes for future debugging or further development                  %
%                                                                         %
% Variables:                                                              %
% r: refers to row index              c: refers to col index              %
% [r/c]old: all matrix indices before the current iteration               %
% [r/c]pre: the matrix indices in the previous iteration                  %
% [r/c]new: the matrix indices in the current iteration                   %
% [r/c]all: the matrix indices containing both [r/c]old and [r/c]new      %
% Index: matrix indices of each iteration (cell of [r/c]all)              %
% MatCal: matrix selected in the current iteration for calculating avg.   %
% MatCalAvg: Mean value of MatCal                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Op = Overlap;
% Deal with strike ranges from 0~45 degrees
if Str >= 0 && Str <= 45
    col = 1; count = 1;
    % How many rows should be added if add one column
    r_inc = 1/tand(Str);
    row = r_inc*col;
    % Set to the size of matrix, if exceeds the matrix width
    if row > size(Grd,1)
        row = size(Grd,1);
    end
    c = col.*ones(1,ceil(row));
    cold{count} = c; cpre = c;
    r = 1:ceil(row);
    rold = r; rpre = r;
    Index{count} = sub2ind(size(Grd),r,c);
    MatCal = Grd(Index{count});
    while length(MatCal) < size(Grd,1)*size(Grd,2)
        count = count + 1;
        col = col + 1;
        if col <= size(Grd,2)
            cnew = cpre + 1;
            row = r_inc*col;
            % Set to the size of matrix, if exceeds the matrix width
            if row > size(Grd,1)
                row = size(Grd,1);
            end
            rnew = 1:ceil(row);
            cadd = length(rnew) - length(rpre);
            rpre = rnew;
            cnew = [cnew,ones(1,cadd)];
            cpre = cnew;
            call = [cold{count-1},cnew];
            rall = [rold,rnew];
            Index{count} = sub2ind(size(Grd),rall,call);
            MatCal = Grd(Index{count});
            MatCalAvg(count) = mean(MatCal,'omitnan');
        elseif col > size(Grd,2)
            outind = cpre == size(Grd,2);
            cpre(outind) = [];
            rpre(outind) = [];
            rnew = rpre;
            cnew = cpre + 1;
            cpre = cnew;
            call = [cold{count-1},cnew];
            rall = [rold,rnew];
            Index{count} = sub2ind(size(Grd),rall,call);
            MatCal = Grd(Index{count});
            MatCalAvg(count) = mean(MatCal,'omitnan');
        end
        cold{count} = call;
        rold = rall;
    end
    % Find the peak of mean value which separates hanging and foot wall
    if PtoN == 1
        Neighbor = 1;
        [~,Mindex] = findlocalmax(MatCalAvg,Neighbor);
        while length(Mindex) > 1
            Neighbor = Neighbor + 1;
            [~,Mindex] = findlocalmax(MatCalAvg,Neighbor);
        end
    elseif PtoN == 0
        Neighbor = 1;
        [~,Mindex] = findlocalmin(MatCalAvg,Neighbor);
        while length(Mindex) > 1
            Neighbor = Neighbor + 1;
            [~,Mindex] = findlocalmin(MatCalAvg,Neighbor);
        end
    end
    LeftIndex = Index{Mindex+Op};
    Left = zeros(size(Grd));
    for i = 1:size(Grd,1)*size(Grd,2)
        if any(i == LeftIndex)
            [r,c] = ind2sub(size(Grd),i);
            Left(r,c) = Grd(r,c);
        else
            [r,c] = ind2sub(size(Grd),i);
            Left(r,c) = nan;
        end
    end
    Edge = setdiff(LeftIndex,Index{Mindex-Op});
    RightIndex = Index{count};
    RightIndex = setdiff(RightIndex,LeftIndex);
    RightIndex = [Edge,RightIndex];
    Right = zeros(size(Grd));
    for i = 1:size(Grd,1)*size(Grd,2)
        if any(i == RightIndex)
            [r,c] = ind2sub(size(Grd),i);
            Right(r,c) = Grd(r,c);
        else
            [r,c] = ind2sub(size(Grd),i);
            Right(r,c) = nan;
        end
    end
    SideA = Left;
    SideB = Right;
    
% Deal with strike ranges from 45~90 degrees
elseif Str > 45 && Str <= 90
    row = 1; count = 1;
    % How many rows should be added if add one column
    c_inc = tand(Str);
    col = c_inc*row;
    % Set to the size of matrix, if exceeds the matrix width
    if col > size(Grd,2)
        col = size(Grd,2);
    end
    r = row.*ones(1,ceil(col));
    rold{count} = r; rpre = r;
    c = 1:ceil(col);
    cold = c; cpre = c;
    Index{count} = sub2ind(size(Grd),r,c);
    MatCal = Grd(Index{count});
    while length(MatCal) < size(Grd,1)*size(Grd,2)
        count = count + 1;
        row = row + 1;
        if row <= size(Grd,1)
            rnew = rpre + 1;
            col = c_inc*row;
            % Set to the size of matrix, if exceeds the matrix width
            if col > size(Grd,2)
                col = size(Grd,2);
            end
            cnew = 1:ceil(col);
            radd = length(cnew) - length(cpre);
            cpre = cnew;
            rnew = [rnew,ones(1,radd)];
            rpre = rnew;
            rall = [rold{count-1},rnew];
            call = [cold,cnew];
            Index{count} = sub2ind(size(Grd),rall,call);
            MatCal = Grd(Index{count});
            MatCalAvg(count) = mean(MatCal,'omitnan');
        elseif row > size(Grd,1)
            outind = rpre == size(Grd,1);
            rpre(outind) = [];
            cpre(outind) = [];
            cnew = cpre;
            rnew = rpre + 1;
            rpre = rnew;
            rall = [rold{count-1},rnew];
            call = [cold,cnew];
            Index{count} = sub2ind(size(Grd),rall,call);
            MatCal = Grd(Index{count});
            MatCalAvg(count) = mean(MatCal,'omitnan');
        end
        rold{count} = rall;
        cold = call;
    end
    % Find the peak of mean value which separates hanging and foot wall
    if PtoN == 1
        Neighbor = 1;
        [~,Mindex] = findlocalmax(MatCalAvg,Neighbor);
        while length(Mindex) > 1
            Neighbor = Neighbor + 1;
            [~,Mindex] = findlocalmax(MatCalAvg,Neighbor);
        end
    elseif PtoN == 0
        Neighbor = 1;
        [~,Mindex] = findlocalmin(MatCalAvg,Neighbor);
        while length(Mindex) > 1
            Neighbor = Neighbor + 1;
            [~,Mindex] = findlocalmin(MatCalAvg,Neighbor);
        end
    end
    LeftIndex = Index{Mindex+Op};
    Left = zeros(size(Grd));
    for i = 1:size(Grd,1)*size(Grd,2)
        if any(i == LeftIndex)
            [r,c] = ind2sub(size(Grd),i);
            Left(r,c) = Grd(r,c);
        else
            [r,c] = ind2sub(size(Grd),i);
            Left(r,c) = nan;
        end
    end
    Edge = setdiff(LeftIndex,Index{Mindex-Op}); 
    RightIndex = Index{count};
    RightIndex = setdiff(RightIndex,LeftIndex);
    RightIndex = [Edge,RightIndex];
    Right = zeros(size(Grd));
    for i = 1:size(Grd,1)*size(Grd,2)
        if any(i == RightIndex)
            [r,c] = ind2sub(size(Grd),i);
            Right(r,c) = Grd(r,c);
        else
            [r,c] = ind2sub(size(Grd),i);
            Right(r,c) = nan;
        end
    end
    SideA = Left;
    SideB = Right;
    SideAIndex = LeftIndex;
    SideBIndex = RightIndex;
end

if fl == 1
    SideA = fliplr(SideA);
    SideB = fliplr(SideB);
    % Flip indices
    Lefttmp = zeros(size(Grd));
    Righttmp = zeros(size(Grd));
    for i = 1:size(Grd,2)
        if mod(size(Grd,2),2) == 1
            Midind = ceil(size(Grd,2)/2);
            Num = i*size(Grd,1)-(size(Grd,1)-1):i*size(Grd,1);
            if i < Midind 
                for j = 1:length(Num)
                    add = (size(Grd,2)-(i*2-1))*size(Grd,1);
                    tmp = Num(j);
                    Ltmp = LeftIndex(tmp == LeftIndex) + add;
                    Rtmp = RightIndex(tmp == RightIndex) + add;
                    if isempty(Ltmp)
                        Ltmp = nan;
                    end
                    if isempty(Rtmp)
                        Rtmp = nan;
                    end
                    Lefttmp(j,i) = Ltmp;
                    Righttmp(j,i) = Rtmp;
                end
            elseif i == Midind
                for j = 1:length(Num)
                    tmp = Num(j);
                    Ltmp = LeftIndex(tmp == LeftIndex);
                    Rtmp = RightIndex(tmp == RightIndex);
                    if isempty(Ltmp)
                        Ltmp = nan;
                    end
                    if isempty(Rtmp)
                        Rtmp = nan;
                    end
                    Lefttmp(j,i) = Ltmp;
                    Righttmp(j,i) = Rtmp;
                end
            elseif i > Midind
                for j = 1:length(Num)
                    tmp = Num(j);
                    sub = (size(Grd,2)-(i*2-1))*size(Grd,1);
                    Ltmp = LeftIndex(tmp == LeftIndex) - sub;
                    Rtmp = RightIndex(tmp == RightIndex) - sub;
                    if isempty(Ltmp)
                        Ltmp = nan;
                    end
                    if isempty(Rtmp)
                        Rtmp = nan;
                    end
                    Lefttmp(j,i) = Ltmp;
                    Righttmp(j,i) = Rtmp;
                end
            end
        elseif mod(size(Grd,2),2) == 0
            Midind = ceil(size(Grd,2)/2);
            Num = i*size(Grd,1)-(size(Grd,1)-1):i*size(Grd,1);
            if i <= Midind
                for j = 1:length(Num)
                    add = (size(Grd,2)-(i*2-1))*size(Grd,1);
                    tmp = Num(j);
                    Ltmp = LeftIndex(tmp == LeftIndex) + add;
                    Rtmp = RightIndex(tmp == RightIndex) + add;
                    if isempty(Ltmp)
                        Ltmp = nan;
                    end
                    if isempty(Rtmp)
                        Rtmp = nan;
                    end
                    Lefttmp(j,i) = Ltmp;
                    Righttmp(j,i) = Rtmp;
                end
            elseif i > Midind
                for j = 1:length(Num)
                    sub = (size(Grd,2)-(i*2-1))*size(Grd,1);
                    tmp = Num(j);
                    Ltmp = LeftIndex(tmp == LeftIndex) - sub;
                    Rtmp = RightIndex(tmp == RightIndex) - sub;
                    if isempty(Ltmp)
                        Ltmp = nan;
                    end
                    if isempty(Rtmp)
                        Rtmp = nan;
                    end
                    Lefttmp(j,i) = Ltmp;
                    Righttmp(j,i) = Rtmp;
                end
            end
        end
    end
    Lefttmp = Lefttmp(:);
    SideAIndex = Lefttmp(~isnan(Lefttmp));
    Righttmp = Righttmp(:);
    SideBIndex = Righttmp(~isnan(Righttmp));
end


end