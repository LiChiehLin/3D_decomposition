%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Lin,Li-Chieh                                %
%                         Department of Geography                         %
%                       National Taiwan University                        %
%                              2021.04.13                                 %
%                                                                         %
% Updates:                                                                %
% (2022.06.15) Add inversion input counts as variable 'Outcount'          %
% (2023.08.21) Add standard error of the inverted displacement for each   %
%              grid                                                       %
%                                                                         %
% Perform 3D displacement inversion on InSAR results                      %
%                                                                         %
% Input:                                                                  %
% 1. InGrd: A cell of Grds (InGrd = {Grd1,Grd2,Grd3...})                  %
% 2. Azimuth: Azimuth angle (Azimuth = [Azi1,Azi2,Azi3...])               %                      
% 3. LookAngle: Incident angle (LookAngle = [Theta1,Theta2,Theta3...])    %
% 4. DispType: LOS or Azi (DispType = {'LOS','Azi','LOS'...})             %
% 5. Orbit: Ascending or Descending (Orbit = {'Asc','Des','Asc'...})      %
% 6. Num: How many inputs (3 or 4 ....)                                   %
%                                                                         %
% Note that the sign of every input should be                             %
% Positive: away from satellite                                           %
% Negative: close to satellite                                            %
%                                                                         %
% Output:                                                                 %
% 1. Out: A cell contains e,n,u disp. in each cell                        %
%    Note that the order is (e,n,u)                                       %
% 2. Outcount: A matrix contains the input counts for the 3D inversion    %
%    The first digit is the azimuth disp. the second is the LOS disp.     %
%    e.g. '21' would be assign to the grid if this grid was inverted from %
%    2 azimuth disp. and 1 LOS disp.                                      %
% 3. OutStderr: A cell contains the standard error of each inverted disp. %
%    for every grid                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Out,OutStderr,Outcount] = InSAR3Ddisp(InGrd,Azimuth,LookAngle,DispType,Orbit,Num)

Row = size(InGrd{1},1);
Col = size(InGrd{1},2);

%Construct G and d matrices
w = ones(Num,1);
G = ones(Num,3);
d = cell(Row,Col);
m = cell(Row,Col);
stderr = cell(Row,Col);
count = zeros(Row,Col);

for i = 1:Num
    %G matrix
    Azi = Azimuth(i);
    Theta = LookAngle(i);
    Type = DispType{i};
    Orb = Orbit{i};
    if strcmp(Type,'Azi')
        if strcmp(Orb,'Asc')
            w(i) = 1;
            G(i,:) = [-sind(Azi),cosd(Azi),0];
        elseif strcmp(Orb,'Des')
            w(i) = 1;
            G(i,:) = [-sind(Azi),-cosd(Azi),0];
        else
            error('Orbit input should be Asc or Des')
        end
    elseif strcmp(Type,'LOS')
        if strcmp(Orb,'Asc')
            w(i) = 1;
            G(i,:) = [cosd(Azi)*sind(Theta),sind(Azi)*sind(Theta),-cosd(90-Theta)];
        elseif strcmp(Orb,'Des')
            w(i) = 1;
            G(i,:) = [-cosd(Azi)*sind(Theta),sind(Azi)*sind(Theta),-cosd(90-Theta)];
        else
            error('Orbit input should be Asc or Des')
        end
    else
        error('Cannot match displacement type: ABORT!')
    end
    
    Wmat = diag(w);
    %d matrix
    Grd = InGrd{i};
    for r = 1:Row
        for c = 1:Col
            tex = ['Construct d Matrix ','Grd:',num2str(i),' row:',num2str(r),' col:',num2str(c)];
            disp(tex)
            d{r,c}(i,1) = Grd(r,c);
        end
    end
end

%Invert to 3D displacement
for r = 1:Row
    for c = 1:Col
        dCal = d{r,c};
        tex = ['Inverting: ','row:',num2str(r),' col:',num2str(c)];
        disp(tex)
        
        Typetmp = DispType(~isnan(dCal));
        Wtmp = Wmat(~isnan(dCal),~isnan(dCal));
        Gtmp = G(~isnan(dCal),:);
        dCaltmp = dCal(~isnan(dCal));
        if (Num - sum(isnan(dCal))) < 2
            m{r,c} = nan;
            count(r,c) = nan;
            stderr{r,c} = nan;
        elseif (Num - sum(isnan(dCal))) == 2
            if strcmp(Typetmp{1},'Azi') && strcmp(Typetmp{2},'Azi')
                Gnew = Wtmp*Gtmp;
                mtmp = Gnew\dCaltmp;
                mtmp(3) = nan;
                m{r,c} = mtmp;
                count(r,c) = 20;
                % Stndard error calculation
                e = Gnew*mtmp - dCaltmp;
                sigma = (e'*e)/length(dCaltmp);
                stderrtmp = sqrt(sigma*diag(inv(Gnew'*Gnew)));
                stderr{r,c} = stderrtmp;
            elseif strcmp(Typetmp{1},'LOS') && strcmp(Typetmp{2},'LOS')
                Gtmp(:,2) = 0;
                Gnew = Gtmp;
                mtmp = Gtmp\dCaltmp;
                mtmp(2) = nan;
                m{r,c} = mtmp;
                count(r,c) = 02;
                % Stndard error calculation
                e = Gnew*mtmp - dCaltmp;
                sigma = (e'*e)/length(dCaltmp);
                stderrtmp = sqrt(sigma*diag(inv(Gnew'*Gnew)));
                stderr{r,c} = stderrtmp;
            else
                m{r,c} = nan;
                count(r,c) = nan;
                stderr{r,c} = nan;
            end
        elseif (Num - sum(isnan(dCal))) > 2
            if sum(strcmp(Typetmp,'LOS')) == length(Typetmp)
                Gtmp(:,2) = 0;
                Gnew = Wtmp*Gtmp;
                mtmp = Gnew\dCaltmp;
                m{r,c} = mtmp;
                count(r,c) = 0 + length(dCaltmp);
                % Stndard error calculation
                e = Gnew*mtmp - dCaltmp;
                sigma = (e'*e)/length(dCaltmp);
                stderrtmp = sqrt(sigma*diag(inv(Gnew'*Gnew)));
                stderr{r,c} = stderrtmp;
            elseif sum(strcmp(Typetmp,'Azi')) == length(Typetmp)
                Gnew = Wtmp*Gtmp;
                mtmp = Gnew\dCaltmp;
                mtmp(3) = nan;
                m{r,c} = mtmp;
                count(r,c) = length(dCaltmp)*10 + 0;
                % Stndard error calculation
                e = Gnew*mtmp - dCaltmp;
                sigma = (e'*e)/length(dCaltmp);
                stderrtmp = sqrt(sigma*diag(inv(Gnew'*Gnew)));
                stderr{r,c} = stderrtmp;
            else
                Gnew = Wtmp*Gtmp;
                mtmp = Gnew\dCaltmp;
                m{r,c} = mtmp;
                count(r,c) = sum(strcmp(Typetmp,'Azi')*10) + sum(strcmp(Typetmp,'LOS'));
                % Stndard error calculation
                e = Gnew*mtmp - dCaltmp;
                sigma = (e'*e)/length(dCaltmp);
                stderrtmp = sqrt(sigma*diag(inv(Gnew'*Gnew)));
                stderr{r,c} = stderrtmp;
            end
        end
    end
end
Out = m;
OutStderr = stderr;
Outcount = count;
end


