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
ind1 = 1;

%plottingAxis = [ -0.3171    1.26   -0.0245    0.5430];
%plottingAxis = [ -0.3171    15.00  -0.0245    2.00];
plottingAxis =[ -0.01 3.5 -0.0245    1.4];
plottingParams.axis = plottingAxis;

ind1 = 190
while ind1<=nFiles
    
    A = dir([pathToDataFileStorage '*.txt']);
    nFiles = length(A);
    
    ind1;
    tStep = ind1*storageStride;
    tTime = tStep*dt;
    
    s = readInDataFile([pathToDataFileStorage A(ind1).name]);
    I = find(s(:,7)==7);
    I1= find(s(:,7)==2);
    Freeparticles = s(I,1:2);
    %AllParticles  = s(:,1:2);
    Velocities    = s(I,3:4);
    Constrained   = s(I1,1:2); 
    x     = Freeparticles(:,1);
    y     = Freeparticles(:,2);
    vx    = Velocities(:,1);
    vy    = Velocities(:,2);
    vMag = sqrt(vx.^2+vy.^2);
    
    xConstr   =  Constrained(:,1);
    yConstr   =  Constrained(:,2);
    
    pressureFree = s(I,6);
    p_ref = 1000 * 9.81 * 0.6;
    pressure = pressureFree / p_ref; 
 
    axis([plottingAxis]);
    
%     if 1*tStep ==2
   
      
        plot(xConstr,yConstr,'k.')
        hold on
        %scatter(x,y,10,vMag,'filled')
        scatter(x,y,5,pressure,'filled')
        %hold on
        %quiver(x,y,vx,vy,2,'k')
        xlabel('distance (m)')
        ylabel('height (m)')
        title(['t = ' sprintf('%f',(ind1-1)*dt*storageStride / t_ref)],'Fontweight','normal');
        niceFig()
        set(gca,'dataaspectratio',[1 1 1]);
        axis([plottingAxis])
        %axis off
        caxis([0 1.0])
        colormap(jet(1024))
        colorbar
%         break
%     end
     drawnow
     match =["dataFile0000","dataFile000","dataFile00","dataFile0"];
     str = erase(A(ind1).name,match);
     cd(pathToImageStorage);
     print(gcf,'-dpng',[str(1:end-4) ,'.png']);
     display([num2str(ind1) ' / ' num2str(nFiles)]);
     ind1 = ind1+1;
end
