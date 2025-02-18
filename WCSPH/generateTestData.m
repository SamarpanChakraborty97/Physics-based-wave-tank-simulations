close all
clear all


%%%%%%%%%%%%%%%%
%RELETIVE PATHS%
%%%%%%%%%%%%%%%%
pathToSPHCode = '../../Release/';
%pathToSPHCode = '../';
pathToHere = mfilename('fullpath');
pathToDataFileStorage = 'dataFiles/';
nameOfInputFile       = 'SPHinputFileInitialization.txt';


%%%%%%%%%%%%%%%%
%FREE PARTICLES%
%%%%%%%%%%%%%%%%
h = 5.0*10^(-2);


%define free particle positions
s.free.xSpan  = [2.0*h 0.25+h];
s.free.ySpan  = [2.0*h 0.3+2.0*h];
xnSpan        = length([2.0*h:1.6*h:0.25+2.0*h])
ynSpan        = length([2.0*h:1.6*h:0.30+2.0*h])
%pause
s.free.nSpan  = [xnSpan ynSpan];


%define free particle characteristice
%if one number, constant
%if three number, mean and standard deviation of gaussian dstribution
%gaussian = 1; 

%parameters from book:
s.free.mass            = 1000*pi*(0.88*h)^2; %mass                         [kg/m^2]
s.free.smoothingLength = h;          %smoothingLength - [m]

s.free.vf          = 1.5;   %maximum fluid velocity
%sprinf('B = ');
s.free.mu          = 0.01; %viscosity
s.free.muBeta      = -2;  %muBeta
s.free.epsilon     = 0.25;  %epsilon in XSPH
s.free.rRef        = 1000; %reference density

s.free.color       = 1;                                      

%%%%%%%%%%%%%%%%%%%%%%%
%CONSTRAINED PARTICLES%
%%%%%%%%%%%%%%%%%%%%%%%

%define container parameters
%container parameters are fixed
s.cons.color = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTATIONAL PARAMETERS%
%%%%%%%%%%%%%%%%%%%%%%%%%%
s.comp.gravity        = -9.8;

%time parameters
s.comp.nT             = 1000*100;
s.comp.dt             = 2*10^(-5);

%export parameters
s.comp.storageStride  = 100;


%%%%%%%%%%%%%%%%
%END USER INPUT%
%%%%%%%%%%%%%%%%


%create paths
I = regexp(pathToHere,'\');
pathToHere = pathToHere(1:I(end));
s.path.pathToCase = ('.\');
s.path.nameOfInputFile  = nameOfInputFile;
pathToCCode       = [pathToSPHCode];

%make the data structures
%s = makeFreeParticles2(s);
s = makeFreeParticles_load(s);
%s = makeContainer2(s);
s = makeContainer3(s);
s = LoadDensity(s); % problem with loading density when reverting back to makeContainer2 


disp([s.freeParticles.nP+ length(s.containerParticles.pos(:,1))])

%2-312
%566- 695

if 0*and(1,length(s.containerParticles.pos(:,1))>853)
    ix = [2:312];
    plot(s.freeParticles.pos(:,1),s.freeParticles.pos(:,2),'b.')
    hold on
    plot(s.containerParticles.pos(:,1),s.containerParticles.pos(:,2),'k.')
    plot(s.containerParticles.pos(ix,1),s.containerParticles.pos(ix,2),'ro')
    axis equal
end

% return

%ouput text file
% outputTextFile(s)
outputTextFile_loadDensity(s)

%execute C code
I = regexp(pathToCCode,'/'); %replace / with \
pathToCCode(I) = '\'; 

I = regexp(pathToDataFileStorage,'/'); %replace / with \
pathToDataFileStorage(I) = '\'; 


dos(['cd ' pathToCCode '&&SPH2DCUDA ' pathToHere nameOfInputFile ' ' pathToHere pathToDataFileStorage ' &'],'-echo');
%dos(['cd ' pathToCCode '&&SPH2DCPPCudaBinning.exe ../inputFile.txt']);
%dos(['cd ' pathToCCode '&&SPH2DCPPCudaBinning.exe ' pathToHere nameOfInputFile]);


