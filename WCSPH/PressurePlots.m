close all
clear all

pathToHere         = mfilename('fullpath');
I = regexp(pathToHere,'\');
pathToHere = pathToHere(1:I(end));
pathToDataFileStorage = [pathToHere 'dataFiles/'];
pathToImageStorage    = [pathToHere 'pressureContour/'];
storageStride  = 400;
dt             = 1*10^(-5);
height  = 0.6;
g = 9.81;
t_ref = sqrt( height / g);
A  = dir([pathToDataFileStorage '*.txt']);
nFiles = length(A);

plottingAxis =[0 5.38 0 2];
plottingParams.axis = plottingAxis;

A = dir([pathToDataFileStorage '*.txt']);

index = [124, 297, 353, 383, 439];
for ind1 = 1:5
tStep = index(ind1)*storageStride;
tTime = tStep*dt;
 s = readInDataFile([pathToDataFileStorage A(index(ind1)).name]);
    I = find(s(:,7)==7);
    I1= find(s(:,7)==2);
    Freeparticles = s(I,1:2); 
    x     = Freeparticles(:,1);
    y     = Freeparticles(:,2);
    XFree = x/height;
    YFree = y/height;
    %AllParticles  = s(:,1:2);
    Velocities    = s(I,3:4);
    Constrained   = s(I1,1:2);
    xConstr   =  Constrained(:,1);
    yConstr   =  Constrained(:,2);
    XCons = xConstr/height;
    YCons = yConstr/height;
    pressureFree = s(I,6);
    p_ref = 1000 * 9.81 * 0.6;
    pressure = pressureFree / p_ref; 
 
    axis([plottingAxis]);
    subplot(5,1,ind1)
     plot(XCons,YCons,'b.','MarkerSize',2)
        hold on
        %scatter(x,y,10,vMag,'filled')
        scatter(XFree,YFree,1,pressure,'filled')
        %niceFig()
        axis([plottingAxis]);
caxis([0 1.0])

        colormap(jet(1024))
        
end

