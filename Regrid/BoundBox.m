%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                         Department of Geography                         %
%                       National Taiwan University                        %
%                              2021.06.01                                 %
%                                                                         %
% Revised: 2021.07.23                                                     %
%          Set additional variable "Datatype" for different data inputs   %
%                                                                         %
%                                                                         %
% Determine the two boundary points for further Regrid.m                  %
%                                                                         %
% Input:                                                                  %
% 1. InGrd: A cell of Grds (InGrd = {Grd1,Grd2,Grd3...})                  %
% 2. Num: How many inputs (3 or 4 ....)                                   %
% 3. Factor: Multiply the increment if desired to make coarser grid cell  %
%    Set 1: Original increment                                            %
%    Set 2: 2 times the original increment. etc.                          %
% 4. Datatype: Whether the input is raw or processed(merged)              %
%              0: Raw (Direct product from GMTSAR)                        %
%                 InGrd is a string cell {'grd1','grd2'...}               %
%              1: Processed (Merged. Store in .mat file)                  %
%                 InGrd is a big cell {{cell1},{cell2}...}                %
%                                                                         %
% Note: This is made to determine the boundary points for the new grid.   %
%                                                                         %
% Output:                                                                 %
% 1. Bound: A 2-by-2 matrix of the Lower-left and Upper-right point.      %
% 2. Interval: Increment of the new grid (Interval = [xinc,yinc])         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Bound,Interval] = BoundBox(InGrd,Num,Factor,Datatype)
%% Check inputs
%Check input is a cell or not
if ~iscell(InGrd)
    error('Input should be a string cell! ABORT!')
end

%Check the number of inputs are consistent
if length(InGrd)~=Num
    error('Num should be consistent with the number of input grds! ABORT!')
end

%% Program start
mul = Factor;
if Datatype == 0
    disp('Inputs are products directly from GMTSAR')
    for i = 1:Num
        file = InGrd{i};
        [x,y,z] = grdread2(file);
        %Increment
        x_inc(i) = x(2) - x(1);
        y_inc(i) = y(2) - y(1);
        %Boundary lon lat
        x_min(i) = x(1);
        x_max(i) = x(end);
        y_min(i) = y(1);
        y_max(i) = y(end);
    end
    Interval = [max(x_inc)*mul;max(y_inc)*mul];

    Bound = [max(x_min),max(y_min);min(x_max),min(y_max)];

elseif Datatype == 1
    disp('Inputs are processed (Merged). Should be stored in .mat file')
    for i = 1:Num
        file = InGrd{i};
        x = file{1};
        y = file{2};
        z = file{3};
        x_inc(i) = x(1,2) - x(1,1);
        y_inc(i) = y(1,1) - y(2,1);
        %Boundary lon lat
        x_min(i) = x(1,1);
        x_max(i) = x(1,end);
        y_min(i) = y(end,1);
        y_max(i) = y(1,1);
    end
    Interval = [max(x_inc)*mul;max(y_inc)*mul];

    Bound = [max(x_min),max(y_min);min(x_max),min(y_max)];
end
end