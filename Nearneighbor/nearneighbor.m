%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                         Department of Geography                         %
%                       National Taiwan University                        %
%                              2022.10.24                                 %
%                                                                         %
%                                                                         %
% Perform 2D interpolation with near-neighbor algorithm                   %
%                                                                         %
%                                                                         %
% Input:                                                                  %
% 1. InMat: Input vector or matrix                                        %
% 2. Radius: Set the radius for searching neighboring values              %
%    e.g. Set 3 will only search up to 3x3 neighboring values             % 
% 3. Weight: Set 1 or 0 to decide whether put weighting or not            %
%    Set 1: Nearer values will have higher weighting (Smoothing)          %
%    Set 0: Take average of all neighboring values                        %
%                                                                         %
% Output:                                                                 %
% 1. OutMat: Interpolated vector or matrix                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OutMat = nearneighbor(InMat,Radius,Weight)
% Find total missing values and indices
missind = find(isnan(InMat(:)));
[missrow,misscol] = ind2sub(size(InMat),missind); % Get row and col index
miss = sum(isnan(InMat(:)));
disp(strcat('Missing values:',num2str(miss)))

OutMat = InMat;
%%%% Search neighbors for each missing values
for i = 1:length(missind)
    fprintf('Working on #%d\n',i);
% Make sure search box is within matrix
    Matind = missind(i);
    row = missrow(i);
    rtmp = (row-Radius):1:(row+Radius);
    col = misscol(i);
    ctmp = (col-Radius):1:(col+Radius);
    rowtmp = rtmp((rtmp>0) & (rtmp <= size(InMat,1))); % rows to be searched
    coltmp = ctmp((ctmp>0) & (ctmp <= size(InMat,2))); % cols to be searched
    
%%%% Start searching neighbors
% Get neighboring values
    rowind = repmat(rowtmp,length(coltmp),1);
    rowind = rowind(:);
    colind = repmat(coltmp,1,length(rowtmp));
    colind = colind(:);
    
    linind = sub2ind(size(InMat),rowind,colind);
    
%%%% No weighting
    if Weight == 0
        neival = InMat(linind);
        OutMat(Matind) = mean(neival,'omitnan');
    elseif Weight == 1
%%%% Include weighting
        Smseq = 1:1:Radius;
        rlow = row - Smseq;
        rlow(rlow <= 0) = 1;
        rup = row + Smseq;
        rup(rup >= max(rowind)) = max(rowind);
        clow = col - Smseq;
        clow(clow <= 0) = 1;
        cup = col + Smseq;
        cup(cup >= max(colind)) = max(colind);
        
        % Extend matrix to the same sizes
        rowindsm = rowind*ones(1,length(rlow));
        rlowsm = ones(length(rowind),1)*rlow;
        rupsm = ones(length(rowind),1)*rup;
        
        colindsm = colind*ones(1,length(clow));
        clowsm = ones(length(colind),1)*clow;
        cupsm = ones(length(colind),1)*cup;
        
        
        % Get indices of neighboring values for each radius
        Smrow = (rowindsm >= rlowsm) & (rowindsm <= rupsm);
        Smcol = (colindsm >= clowsm) & (colindsm <= cupsm);
        Smlinind = linind*ones(1,Radius);
        Smind = Smlinind(Smrow & Smcol);
        Smneival = InMat(Smind);
        
        % Put value to the interpolated matrix
        OutMat(Matind) = mean(Smneival,'omitnan');
    else
        error('Weight should be 1 or 0 (Yes or No)');
    end
    
end
 
end