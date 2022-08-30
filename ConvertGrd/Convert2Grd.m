%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                         Department of Geography                         %
%                       National Taiwan University                        %
%                              2021.07.23                                 %
%                                                                         %
% Ver2. 2021.08.30                                                        %
% Make it applicable for converting only one grd                          %
%                                                                         %
%                                                                         %
% This is to convert 3D displacement inversion to grd and txt             %
%                                                                         %
% Input:                                                                  %
% 1. Inv: A cell contain five columns (Lon,Lat,E,N,U)                     %
%         All five columns have the same size                             %
% 2. Inc: Increment of X and Y directions                                 %
%         Can take the result from BoundBox.m                             %
% 3. Step: Whether is one-step or two-step. A string (step = 'one')       %
%          If converting only one grd, then Step is the output grd name   %
%                                                                         %
% Note: You still have to go to terminal and csh the .csh created by this %
%       code!!                                                            %
%                                                                         %
% Output:                                                                 %
% 1. The command you should type in the terminal                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = Convert2Grd(Inv,Inc,Step)
if length(Inv) == 3;
disp('Input cell should be 3*1 Lon Lat Z')

%
% Set the output file name
%
OutName = Step;
% 
% Get each column
%
Lon = Inv{1};
Lat = Inv{2};
Z = Inv{3};
%
% Output
%
Zout = [Lon(:),Lat(:),Z(:)];

fid = fopen(strcat(OutName,'.txt'),'w');
[m,n] = size(Zout);
for i = 1:m
    for j = 1:n
        fprintf(fid,'%f %f %f',Zout(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);


%%%%%%%%%%%%%%%
% GMT xyz2grd %
%%%%%%%%%%%%%%%
Rxmin = Lon(1,1);
Rxmax = Lon(1,end);
Rymin = Lat(end,1);
Rymax = Lat(1,1);

Ix = Inc(1);
Iy = Inc(2);

% 32 represent "space"
commandZ = strcat('echo',32,'gmt',32,'xyz2grd',32,OutName,'.txt',32,'-G',OutName,'.grd',32,...
    '-R',num2str(Rxmin),'/',num2str(Rxmax),'/',num2str(Rymin),'/',...
    num2str(Rymax),32,'-I',num2str(Ix),'/',num2str(Iy),32,'-V',32,'>',32,'XYZ2GRD.csh');


system(commandZ);
xyz2grd = 'csh XYZ2GRD.csh';
result = xyz2grd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif length(Inv) == 5
    disp('Input cell should be 5*1 Lon Lat E N U')

%
% Set the output file name
%
Ename = strcat('E_',Step);
Nname = strcat('N_',Step);
Uname = strcat('U_',Step);
%
% Get each column
%
Lon = Inv{1};
Lat = Inv{2};
E = Inv{3};
N = Inv{4};
U = Inv{5};

%
% Output
%
Eout = [Lon(:),Lat(:),E(:)];
Nout = [Lon(:),Lat(:),N(:)];
Uout = [Lon(:),Lat(:),U(:)];

fid = fopen(strcat(Ename,'.txt'),'w');
[m,n] = size(Eout);
for i = 1:m
    for j = 1:n
        fprintf(fid,'%f %f %f',Eout(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);


fid = fopen(strcat(Nname,'.txt'),'w');
[m,n] = size(Nout);
for i = 1:m
    for j = 1:n
        fprintf(fid,'%f %f %f',Nout(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);


fid = fopen(strcat(Uname,'.txt'),'w');
[m,n] = size(Uout);
for i = 1:m
    for j = 1:n
        fprintf(fid,'%f %f %f',Uout(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

%%%%%%%%%%%%%%%
% GMT xyz2grd %
%%%%%%%%%%%%%%%
Rxmin = Lon(1,1);
Rxmax = Lon(1,end);
Rymin = Lat(end,1);
Rymax = Lat(1,1);

Ix = Inc(1);
Iy = Inc(2);

% 32 represent "space"
commandE = strcat('echo',32,'gmt',32,'xyz2grd',32,Ename,'.txt',32,'-G',Ename,'.grd',32,...
    '-R',num2str(Rxmin),'/',num2str(Rxmax),'/',num2str(Rymin),'/',...
    num2str(Rymax),32,'-I',num2str(Ix),'/',num2str(Iy),32,'-V','>',32,'XYZ2GRD.csh');

commandN = strcat('echo',32,'gmt',32,'xyz2grd',32,Nname,'.txt',32,'-G',Nname,'.grd',32,...
    '-R',num2str(Rxmin),'/',num2str(Rxmax),'/',num2str(Rymin),'/',...
    num2str(Rymax),32,'-I',num2str(Ix),'/',num2str(Iy),32,'-V','>>',32,'XYZ2GRD.csh');

commandU = strcat('echo',32,'gmt',32,'xyz2grd',32,Uname,'.txt',32,'-G',Uname,'.grd',32,...
    '-R',num2str(Rxmin),'/',num2str(Rxmax),'/',num2str(Rymin),'/',...
    num2str(Rymax),32,'-I',num2str(Ix),'/',num2str(Iy),32,'-V','>>',32,'XYZ2GRD.csh');

system(commandE);
system(commandN);
system(commandU);
xyz2grd = 'csh XYZ2GRD.csh';
result = xyz2grd;
end
end


