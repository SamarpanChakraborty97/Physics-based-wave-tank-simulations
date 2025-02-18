close all
clear all

pathToSPHCode = 'D:/WCSPH_3D/WCSPH_3D/x64/Release/Dam_Break/';
%pathToSPHCode = '../';
pathToHere = mfilename('fullpath');
pathToDataFileStorage = 'dataFiles/';
nameOfInputFile       = 'SPHinputFileInitialization.txt';


%%%%%%%%%%%%%%%%
%FREE PARTICLES%
%%%%%%%%%%%%%%%%
height = 3 * 10^(-1); % height of water column
width = 6.1 * 10^(-1); % width of water column
length1 = 4 * 10^(-1); % length of water coulmn
length2 = 12 * 10^(-1);% length of second water basin

dx = height/10; % particle spacing = dx;
h = 2 * dx;

s.free.xSpan1=[h length1-h];
s.free.ySpan1=[h width-h];
s.free.zSpan1=[h height];

s.free.xSpan2=[length1+h 0.9-dx];
s.free.ySpan2 =[h width-h];
s.free.zSpan2=[h 0.01];

s.free.xSpan3=[0.9 1.02];
s.free.ySpan3=[h 0.25-dx];
s.free.zSpan3 =[h 0.01];

s.free.xSpan4 = [0.9 1.02];
s.free.ySpan4=[0.37+dx width-h];
s.free.zSpan4 =[h 0.01];

s.free.xSpan5 = [1.02+dx 1.6-h];
s.free.ySpan5 = [h width-h];
s.free.zSpan5 = [h 0.01];

xnSpan1=length(h:dx:length1-h);
ynSpan1=length(h:dx:width-h);
znSpan1=length(h:dx:height);

xnSpan2=length(length1+h:dx:0.9-dx);
ynSpan2 = ynSpan1;
znSpan2=length(h:dx:0.01);

xnSpan3=length(0.9:dx:1.02);
ynSpan3 = length(h:dx:0.25-dx);
znSpan3=znSpan2;

xnSpan4=xnSpan3;
ynSpan4 = length(0.37+dx:dx:width-h);
znSpan4=znSpan2;

xnSpan5=length(1.02+dx:dx:1.6-h);
ynSpan5=ynSpan1;
znSpan5=znSpan2;

s.free.nSpan1 = [xnSpan1 ynSpan1 znSpan1];
s.free.nSpan2 = [xnSpan2 ynSpan2 znSpan2];
s.free.nSpan3 = [xnSpan3 ynSpan3 znSpan3];
s.free.nSpan4 = [xnSpan4 ynSpan4 znSpan4];
s.free.nSpan5 = [xnSpan5 ynSpan5 znSpan5];

rho_a = 1000;
Va = dx ^ 3;
m = rho_a * Va;
s.free.mass            = m;       %mass  [kg/m^3] 
s.free.smoothingLength = h;          %smoothingLength - [m]
s.free.vf          = 3.0;   %maximum fluid velocity
s.free.nu          = 1 * 10^(-6);
s.free.delta       = 0.1; %kinematic viscosity
s.free.muAlpha     = 0.1; %-0.4
s.free.muBeta      = 5.0;  %muBeta -0.5
s.free.epsilon     = 0;  %epsilon in XSPH %0.25
s.free.rRef        = 1000; %reference density
s.free.color       = 1;

%%%%%%%%%%%%%%%%%%%%%%%
%CONSTRAINED PARTICLES%
%%%%%%%%%%%%%%%%%%%%%%%

%define container parameters
%container parameters are fixed
s.cons.color = 2;
%s.measured.color = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTATIONAL PARAMETERS%
%%%%%%%%%%%%%%%%%%%%%%%%%%
s.comp.gravity        = -9.8;

%time parameters
s.comp.nT             = 800 * 300;
s.comp.dt             = 5*10^(-6);

%export parameters
s.comp.storageStride  = 800;

%%%%%%%%%%%%%%%%
%END USER INPUT%
%%%%%%%%%%%%%%%%

I = regexp(pathToHere,'\');
pathToHere = pathToHere(1:I(end));
s.path.pathToCase = ('.\');
s.path.nameOfInputFile  = nameOfInputFile;
pathToCCode       = [pathToSPHCode];

%make the data structures
s = makeFreeParticles2(s);
%s = makeFreeParticles_load(s);
%s = pressurePoint(s);
s = makeContainer3 (s);
%s = makeContainer3_load(s);

plot3(s.freeParticles.pos(:,1),s.freeParticles.pos(:,2),s.freeParticles.pos(:,3),'b.')
hold on
%plot3(s.containerParticles.pos(:,1),s.containerParticles.pos(:,2),s.containerParticles.pos(:,3),'k.')
plot3(s.con.pos(:,1),s.con.pos(:,2),s.con.pos(:,3),'k.')
%hold on;
plot3(s.gate.pos(:,1),s.gate.pos(:,2),s.gate.pos(:,3),'g.')
%hold on;
plot3(s.pillar.pos(:,1),s.pillar.pos(:,2),s.pillar.pos(:,3),'r.')

axis equal
f = gcf;
%saveas(f,'DamBreak2.png')

outputTextFile(s)

disp(s.freeParticles.nP)
disp([s.freeParticles.nP+ length(s.containerParticles.pos(:,1))])

%execute C code
% I = regexp(pathToCCode,'/'); %replace / with \
% pathToCCode(I) = '\'; 
% 
% I = regexp(pathToDataFileStorage,'/'); %replace / with \
% pathToDataFileStorage(I) = '\';

dos(['cd ' pathToCCode '&&WCSPH_3D.exe ' pathToHere nameOfInputFile ' ' pathToHere pathToDataFileStorage ' &'],'-echo');
