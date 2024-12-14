%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Lin,Li-Chieh                                %
%                     Earth and Planetary Sciences                        %
%                  University of California, Riverside                    %
%                              2021.04.13                                 %
%                                                                         %
% Updates:                                                                %
% (2022.06.15) Add inversion input counts as variable 'Outcount'          %
% (2023.08.21) Add standard error of the inverted displacement for each   %
%              grid                                                       %
% (2023.10.25) Make pixel-wise azimuth and incidence angle input          %
% (2024.07.01) MAJOR UPDATES!                                             %
% 1. Replace input "Orbit" with "LookDirection" to accomodate various     %
%    flying direction and look direction                                  %
% 2. Change Azimuth to be w.r.t. north. North is zero and positive along  %
%    clock-wise direction.                                                %
% 3. Remove input "Num"                                                   %
% (2024.07.11) Support input for weighting of each input "Weights"        %
% Three ways for input argument "Weights": (Details see below)            %
% 1. Have the same size of InMat, each pixel has different weights        %
% 2. Put in the same weights for each input for all pixel                 %
% 3. Simply put a string 'no' to assume uniform weighting (No weights)    %
%                                                                         %
% Perform 3D displacement inversion on InSAR results                      %
%                                                                         %
% Input:                                                                  %
% 1. InGrd: A cell of Grds (InGrd = {Grd1,Grd2,Grd3...})                  %
% 2. Weights: Three formats                                               %
%  2.1. A cell of Grds (Same size as "InGrd"). Different weights for each %
% pixel and input                                                         %
%  2.2. A Mx1 or 1xM vector assigning the same weights for the same input.%
% Order should be the same as the inputs.                                 %
%  2.3. A string 'no' to assume uniform weighting (No weights)            %
% 3. Azimuth: Azimuth angle (Azimuth = [Azi1,Azi2,Azi3...]) or            %
%    (Azimuth = {Azi1,Azi2,Azi3}), each is a matrix with same size        %
% 4. LookAngle: Incident angle (LookAngle = [Theta1,Theta2,Theta3...])    %
%    (LookAngle = {Theta1,Theta2,Theta3}), each is a matrix with same size%
% 5. DispType: LOS or Azi (DispType = {'LOS','Azi','LOS'...})             %
% 6. LookDirection: Left or Right (LookDirection = {'l','r','r'...})      %
%                                                                         %
% Note that the sign of every input should be                             %
% Positive: away from satellite                                           %
% Negative: close to satellite                                            %
%                                                                         %
% Output:                                                                 %
% 1. Out: A cell contains e,n,u disp. in each cell                        %
%    Note that the order is (e,n,u)                                       %
% 2. OutModelVar: A cell contains the model variance of each inverted disp%
%    for every grid.                                                      %
% 3. Outcount: A matrix contains the input counts for the 3D inversion    %
%    The first digit is the azimuth disp. the second is the LOS disp.     %
%    e.g. '21' would be assign to the grid if this grid was inverted from %
%    2 azimuth disp. and 1 LOS disp.                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Out,OutModelVar,Outcount] = InSAR3Ddisp(InGrd,Weights,Azimuth,LookAngle,DispType,LookDirection)

Row = size(InGrd{1},1);
Col = size(InGrd{1},2);
Num = length(InGrd);

if iscell(Weights) && (length(Weights) == length(InGrd))
    Wtype = '1';
elseif ~iscell(Weights) && (length(Weights) == length(InGrd))
    Wtype = '2';
elseif strcmp(Weights,'no')
    Wtype = 'no';
else
    error( 'u:stuffed:it' , ['Weight input should be:\n', ...
        '(1) a cell containing the same size of InGrd (Different weights for each pixel and each input), or \n', ...
        '(2) a vector of the weights for each input observation (Different weights for different input same for all pixels), or \n', ...
        '(3) simply put a string "no" to assume uniform weighting']);
end

%Construct G and d matrices
W = cell(Row,Col);
G = cell(Row,Col);
d = cell(Row,Col);
m = cell(Row,Col);
stderr = cell(Row,Col);
count = zeros(Row,Col);
for i = 1:Num
    if iscell(Azimuth) && iscell(LookAngle)
        Azi = Azimuth{i};
        Theta = LookAngle{i};
        Azi_check = mean(Azi(:),'omitnan');
    elseif iscell(Azimuth) && ~iscell(LookAngle)
        Azi = Azimuth{i};
        Theta = LookAngle(i)*ones(Row,Col);
        Azi_check = mean(Azi(:),'omitnan');
    elseif ~iscell(Azimuth) && iscell(LookAngle)
        Azi = Azimuth(i)*ones(Row,Col);
        Theta = LookAngle{i};
        Azi_check = Azi;
    elseif ~iscell(Azimuth) && ~iscell(LookAngle)
        Azi = Azimuth(i)*ones(Row,Col);
        Theta = LookAngle(i)*ones(Row,Col);
        Azi_check = Azi;
    end

    % Which quadrant is the sensor flying
    if Azi_check > 0 && Azi_check < 90
        Quad = '1';
    elseif Azi_check > 90 && Azi_check < 180
        Quad = '4';
    elseif Azi_check > 180 && Azi_check < 270
        Quad = '3';
    elseif Azi_check > 270 && Azi_check < 360
        Quad = '2';
    end

    % Constrcut G and d matrix
    Type = DispType{i};
    LR = LookDirection{i};
    Grd = InGrd{i};
    for r = 1:Row
        for c = 1:Col
            tex = ['Construct G and d Matrix ','Grd:',num2str(i),' row:',num2str(r),' col:',num2str(c)];
            disp(tex)
            alph = Azi(r,c);
            thet = Theta(r,c);
            if strcmp(LR,'l') && strcmp(Quad,'1')
                GLos = [-cosd(alph).*sind(thet),sind(alph).*sind(thet),-cosd(thet)];
                GAzi = [sind(alph),cosd(alph),0];
            elseif strcmp(LR,'r') && strcmp(Quad,'1')
                GLos = [cosd(alph).*sind(thet),-sind(alph).*sind(thet),-cosd(thet)];
                GAzi = [sind(alph),cosd(alph),0];
            elseif strcmp(LR,'l') && strcmp(Quad,'2')
                GLos = [-cosd(alph).*sind(thet),sind(alph).*sind(thet),-cosd(thet)];
                GAzi = [sind(alph),cosd(alph),0];
            elseif strcmp(LR,'r') && strcmp(Quad,'2')
                GLos = [cosd(alph).*sind(thet),-sind(alph).*sind(thet),-cosd(thet)];
                GAzi = [sind(alph),cosd(alph),0];
            elseif strcmp(LR,'l') && strcmp(Quad,'3')
                GLos = [-cosd(alph).*sind(thet),sind(alph).*sind(thet),-cosd(thet)];
                GAzi = [sind(alph),cosd(alph),0];
            elseif strcmp(LR,'r') && strcmp(Quad,'3')
                GLos = [cosd(alph).*sind(thet),-sind(alph).*sind(thet),-cosd(thet)];
                GAzi = [sind(alph),cosd(alph),0];
            elseif strcmp(LR,'l') && strcmp(Quad,'4')
                GLos = [-cosd(alph).*sind(thet),sind(alph).*sind(thet),-cosd(thet)];
                GAzi = [sind(alph),cosd(alph),0];
            elseif strcmp(LR,'r') && strcmp(Quad,'4')
                GLos = [cosd(alph).*sind(thet),-sind(alph).*sind(thet),-cosd(thet)];
                GAzi = [sind(alph),cosd(alph),0];
            else 
                error('LookDirection: l or r')
            end

            % Make final Green's function
            if strcmp(Type,'Azi')
                G{r,c}(i,:) = GAzi;
            elseif strcmp(Type,'LOS')
                G{r,c}(i,:) = GLos;
            else
                error('Cannot match displacement type: ABORT!');
            end

            % Make weighting matrix
            if strcmp(Wtype,'1')
                WGrd = Weights{i};
                W{r,c}(i,i) = WGrd(r,c);
            elseif strcmp(Wtype,'2')
                WGrd = Weights(i);
                W{r,c}(i,i) = WGrd;                
            elseif strcmp(Wtype,'3')
                WGrd = 1;
                W{r,c}(i,i) = WGrd;
            end

            % Make data vector
            d{r,c}(i,1) = Grd(r,c);
        end
    end
end

%Invert to 3D displacement
for r = 1:Row
    for c = 1:Col
        dCal = d{r,c};
        GCal = G{r,c};
        WCal = W{r,c};
        tex = ['Inverting: ','row:',num2str(r),' col:',num2str(c)];
        disp(tex)
        
        Typetmp = DispType(~isnan(dCal));
        Wtmp = WCal(~isnan(dCal),~isnan(dCal));
        Gtmp = GCal(~isnan(dCal),:);
        dCaltmp = dCal(~isnan(dCal));
        if (Num - sum(isnan(dCal))) < 2
            m{r,c} = nan;
            count(r,c) = nan;
            stderr{r,c} = nan;
        elseif (Num - sum(isnan(dCal))) == 2
            if strcmp(Typetmp{1},'Azi') && strcmp(Typetmp{2},'Azi')
                right = transpose(Gtmp)*inv(Wtmp)*Gtmp;
                mtmp = inv(right)*transpose(Gtmp)*inv(Wtmp)*dCaltmp;
                mtmp(3) = nan;
                m{r,c} = mtmp;
                count(r,c) = 20;
                stderr{r,c} = diag(inv(right));
            elseif strcmp(Typetmp{1},'LOS') && strcmp(Typetmp{2},'LOS')
                %Gtmp(:,2) = 0;
                mtmp = Gtmp\dCaltmp;
                mtmp(2) = nan;
                m{r,c} = mtmp;
                count(r,c) = 02;
                stderr{r,c} = diag(inv(transpose(Gtmp)*Gtmp));
            else
                m{r,c} = nan;
                count(r,c) = nan;
                stderr{r,c} = nan;
            end
        elseif (Num - sum(isnan(dCal))) > 2
            if sum(strcmp(Typetmp,'LOS')) == length(Typetmp)
                right = transpose(Gtmp)*inv(Wtmp)*Gtmp;
                mtmp = inv(right)*transpose(Gtmp)*inv(Wtmp)*dCaltmp;
                m{r,c} = mtmp;
                count(r,c) = 0 + length(dCaltmp);
                stderr{r,c} = diag(inv(right));
            elseif sum(strcmp(Typetmp,'Azi')) == length(Typetmp)
                right = transpose(Gtmp)*inv(Wtmp)*Gtmp;
                mtmp = inv(right)*transpose(Gtmp)*inv(Wtmp)*dCaltmp;
                mtmp(3) = nan;
                m{r,c} = mtmp;
                count(r,c) = length(dCaltmp)*10 + 0;
                stderr{r,c} = diag(inv(right));
            else
                right = transpose(Gtmp)*inv(Wtmp)*Gtmp;
                mtmp = inv(right)*transpose(Gtmp)*inv(Wtmp)*dCaltmp;
                m{r,c} = mtmp;
                count(r,c) = sum(strcmp(Typetmp,'Azi')*10) + sum(strcmp(Typetmp,'LOS'));
                stderr{r,c} = diag(inv(right));
            end
        end
    end
end
Out = m;
OutModelVar = stderr;
Outcount = count;
end


