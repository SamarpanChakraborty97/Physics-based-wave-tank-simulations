close all
clear all

pathToHere         = mfilename('fullpath');
I = regexp(pathToHere,'\');
pathToHere = pathToHere(1:I(end));
pathToDataFileStorage = [pathToHere 'dataFiles/'];

storageStride  = 400;
dt             = 1*10^(-5);
height  = 0.6;
g = 9.81;
t_ref = sqrt( height / g);
A  = dir([pathToDataFileStorage '*.txt']);
nFiles = length(A);
ind1 = 1;

ind1 = 1;
P = [];
T=[];
if 1
while ind1<=nFiles
    
    A = dir([pathToDataFileStorage '*.txt']);
    nFiles = length(A);
    
    ind1
    tStep = ind1*storageStride;
    tTime = tStep*dt;
    
    s = readInDataFile([pathToDataFileStorage A(ind1).name]);
    I = find(s(:,7)==5);
%     x     = Freeparticles(:,1%     I1= find(s(:,7)==2);
%     Freeparticles = s(I,1:2);
%     %AllParticles  = s(:,1:2);
%     Velocities    = s(I,3:4);
%     Constrained   = s(I1,1:2); );
%     y     = Freeparticles(:,2);
%     vx    = Velocities(:,1);
%     vy    = Velocities(:,2);
%     vMag = sqrt(vx.^2+vy.^2);
    
%     xConstr   =  Constrained(:,1);
%     yConstr   =  Constrained(:,2);
    
    pressureFree = s(I,6);
    p_ref = 1000 * 9.81 * 0.6;
    pressure = pressureFree / p_ref;
    
    time = tTime / t_ref;
    P(ind1) = pressure;
    T(ind1) = time;
    ind1 = ind1 + 1;
end
end
dataset = xlsread('dataset.xlsx','Default Dataset (1)','A1:B91');
x_valid = dataset(:,1);
y_valid = dataset(:,2);

plot(T,P);
hold on;
plot(x_valid, y_valid, '^')
legend('SPH_2D: H/\Deltax = 50','Experiment')
xlabel('t(g/H)^{1/2}')
ylabel('P/\rho gh')

    
    
    