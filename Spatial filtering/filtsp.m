%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               Lin,Li-Chieh                              %
%                         Department of Geography                         %
%                       National Taiwan University                        %
%                                                                         %
% Ver.1: 2021.04.09                                                       %
% Ver.2: 2021.04.12 Add several more filters and change                   %
% function name to filtsp.m                                               %
%                                                                         %
% Perform spatial smoothing on InSAR disp.(Moving window)                 %
%                                                                         %
% Input:                                                                  %
% 1. Ingrd: Grd that is needed to be smoothed(N by M matrix)              %
% 2. ws: Window size (3, 5, 7 etc.)                                       %
% 3. type:                                                                %
% 'lowp'=low-pass filter                                                  %
% 'highp'=high-pass filter                                                %
% 'gauss'=gaussian filter                                                 %
% 4. sigma: For Gaussian filter(1, 2, 3 etc.)                             %
% 5. nanskip: Ignore nan or not. 1=gives nan if has one in window         %
%                                                                         %
% Example: Out = filtsp(Ingrd,5,'lowp',2)                                 %
%                                                                         %
% Output:                                                                 %
% 1. Outgrd: Grd that is filtered                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Outgrd] = filtsp(Ingrd,ws,type,sigma,nanskip)
if ~ischar(type)
    error('MyComponent:incorrectType',...
        'Error. \n"type" has to be a string variable \n Only accept: "lowp","highp","gauss"');
end

Grd = Ingrd;
%% Construct window according to 'type'
filt = type;
if strcmp(filt,'lowp')
    Win = (1/ws^2)*ones(ws,ws);
elseif strcmp(filt,'highp')
    Win = -1*ones(ws,ws);
    Win(ceil(ws/2),ceil(ws/2)) = (ws^2-1);
    Win = 1/ws^2.*Win;
elseif strcmp(filt,'gauss')
    sigf = sigma;
    GaussMat = -floor(ws/2):floor(ws/2);
    Win = ones(ws,ws);
    row = 0;
    col = 0;
    for i = GaussMat
        row = row + 1;
        for j = GaussMat
            col = col + 1;
            Win(row,col) = (1/(2*pi*sigf^2))*exp(-1*((i^2+j^2)/(2*sigf^2)));
        end
        col = 0;
    end
    Win = Win/sum(sum(Win));
end

%% Boundary grids according to window size
Bound = floor(ws/2);
Outgrd = ones(size(Grd));
naskip = nanskip;
if naskip == 1
    for i = 1:size(Grd,1)
    for j = 1:size(Grd,2)
        %Calculating on uppermost grids
        if (i - Bound <= 0)
            Rowind = unique([ceil(i/Bound):i,i:i+Bound]);
            %Upper-left
            if (j - Bound <= 0)
                Colind = unique([ceil(j/Bound):j,j:j+Bound]);
                MatCal = Grd(Rowind,Colind);
                WinB = Win(Rowind+(ws-max(Rowind)),Colind+(ws-max(Colind)));
                %WinB = Win(1:size(MatCal,1),1:size(MatCal,2));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(sum(MatMul));
            %Upper-right
            elseif (j + Bound > size(Grd,2))
                Colind = unique([j-Bound:j,j:size(Grd,2)]);
                MatCal = Grd(Rowind,Colind);
                WinB = Win(Rowind+(ws-max(Rowind)),1:length(Colind));
                %WinB = Win(1:size(MatCal,1),1:size(MatCal,2));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(sum(MatMul));
            %Elsewhere
            else
                Colind = j-Bound:j+Bound;
                MatCal = Grd(Rowind,Colind);
                WinB = Win(Rowind+(ws-max(Rowind)),1:size(MatCal,2));
                %WinB = Win(1:size(MatCal,1),1:size(MatCal,2));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(sum(MatMul));
            end
            
        %Calculating on lowermost grids    
        elseif (i + Bound >= size(Grd,1)) 
            Rowind = i-Bound:size(Grd,1);
            %Lower-left
            if (j - Bound <= 0)
                Colind = unique([ceil(j/Bound):j,j:j+Bound]);
                MatCal = Grd(Rowind,Colind);
                WinB = Win(1:length(Rowind),Colind+(ws-max(Colind)));
                %WinB = Win(1:size(MatCal,1),1:size(MatCal,2));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(sum(MatMul));
            %Lower-right
            elseif (j + Bound > size(Grd,2))
                Colind = unique([j-Bound:j,j:size(Grd,2)]);
                MatCal = Grd(Rowind,Colind);
                WinB = Win(1:length(Rowind),1:length(Colind));
                %WinB = Win(1:size(MatCal,1),1:size(MatCal,2));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(sum(MatMul));
            %Elsewhere
            else
                Rowind = i-Bound:size(Grd,1);
                Colind = j-Bound:j+Bound;
                MatCal = Grd(Rowind,Colind);
                WinB = Win(1:length(Rowind),1:size(MatCal,2));
                %WinB = Win(1:size(MatCal,1),1:size(MatCal,2));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(sum(MatMul));
            end
        else
            Rowind = i-Bound:i+Bound;
        %Calculating on left-side grids (Skip uppermost and lowermost)
            if (j - Bound <= 0)
                Colind = unique([ceil(j/Bound):j,j:j+Bound]);
                MatCal = Grd(Rowind,Colind);
                WinB = Win(1:size(MatCal,1),Colind+(ws-max(Colind)));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(sum(MatMul));
        %Calculating on right-side grids (Skip uppermost and lowermost)
            elseif (j + Bound > size(Grd,2))
                Colind = unique([j-Bound:j,j:size(Grd,2)]);
                MatCal = Grd(Rowind,Colind);
                WinB = Win(1:size(MatCal,1),1:length(Colind));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(sum(MatMul));
        %Elsewhere of the whole grid
            else
                MatCal = Grd((i-Bound:i+Bound),(j-Bound:j+Bound));
                WinB = Win(1:size(MatCal,1),1:size(MatCal,2));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(sum(MatMul));
            end
        end
    end
    end

elseif naskip == 0
    for i = 1:size(Grd,1)
    for j = 1:size(Grd,2)
        %Calculating on uppermost grids
        if (i - Bound <= 0)
            Rowind = unique([ceil(i/Bound):i,i:i+Bound]);
            %Upper-left
            if (j - Bound <= 0)
                Colind = unique([ceil(j/Bound):j,j:j+Bound]);
                MatCal = Grd(Rowind,Colind);
                WinB = Win(Rowind+(ws-max(Rowind)),Colind+(ws-max(Colind)));
                %WinB = Win(1:size(MatCal,1),1:size(MatCal,2));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(MatMul(~isnan(MatMul)));
            %Upper-right
            elseif (j + Bound > size(Grd,2))
                Colind = unique([j-Bound:j,j:size(Grd,2)]);
                MatCal = Grd(Rowind,Colind);
                WinB = Win(Rowind+(ws-max(Rowind)),1:length(Colind));
                %WinB = Win(1:size(MatCal,1),1:size(MatCal,2));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(MatMul(~isnan(MatMul)));
            %Elsewhere
            else
                Colind = j-Bound:j+Bound;
                MatCal = Grd(Rowind,Colind);
                WinB = Win(Rowind+(ws-max(Rowind)),1:size(MatCal,2));
                %WinB = Win(1:size(MatCal,1),1:size(MatCal,2));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(MatMul(~isnan(MatMul)));   
            end
            
        %Calculating on lowermost grids    
        elseif (i + Bound >= size(Grd,1)) 
            Rowind = i-Bound:size(Grd,1);
            %Lower-left
            if (j - Bound <= 0)
                Colind = unique([ceil(j/Bound):j,j:j+Bound]);
                MatCal = Grd(Rowind,Colind);
                WinB = Win(1:length(Rowind),Colind+(ws-max(Colind)));
                %WinB = Win(1:size(MatCal,1),1:size(MatCal,2));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(MatMul(~isnan(MatMul)));
            %Lower-right
            elseif (j + Bound > size(Grd,2))
                Colind = unique([j-Bound:j,j:size(Grd,2)]);
                MatCal = Grd(Rowind,Colind);
                WinB = Win(1:length(Rowind),1:length(Colind));
                %WinB = Win(1:size(MatCal,1),1:size(MatCal,2));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(MatMul(~isnan(MatMul)));
            %Elsewhere
            else
                Rowind = i-Bound:size(Grd,1);
                Colind = j-Bound:j+Bound;
                MatCal = Grd(Rowind,Colind);
                WinB = Win(1:length(Rowind),1:size(MatCal,2));
                %WinB = Win(1:size(MatCal,1),1:size(MatCal,2));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(MatMul(~isnan(MatMul)));
            end
        else
            Rowind = i-Bound:i+Bound;
        %Calculating on left-side grids (Skip uppermost and lowermost)
            if (j - Bound <= 0)
                Colind = unique([ceil(j/Bound):j,j:j+Bound]);
                MatCal = Grd(Rowind,Colind);
                WinB = Win(1:size(MatCal,1),Colind+(ws-max(Colind)));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(MatMul(~isnan(MatMul)));
        %Calculating on right-side grids (Skip uppermost and lowermost)
            elseif (j + Bound > size(Grd,2))
                Colind = unique([j-Bound:j,j:size(Grd,2)]);
                MatCal = Grd(Rowind,Colind);
                WinB = Win(1:size(MatCal,1),1:length(Colind));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(MatMul(~isnan(MatMul)));
        %Elsewhere of the whole grid
            else
                MatCal = Grd((i-Bound:i+Bound),(j-Bound:j+Bound));
                WinB = Win(1:size(MatCal,1),1:size(MatCal,2));
                MatMul = MatCal.*WinB;
                Outgrd(i,j) = sum(MatMul(~isnan(MatMul)));
            end
        end
    end
    end

else
    error('nanskip should be 1 or 0')
end

end

