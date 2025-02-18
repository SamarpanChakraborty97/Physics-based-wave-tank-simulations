close all
clear all


%x y vx vy xi
% pathToSPHCode         = '\\Client\H$\ENME489ISpring2012\codes\SPH2DCPPCuda\';
% pathToDataFileStorage = '\\Client\H$\ENME489ISpring2012\codes\SPH2DCPPCuda\standardWaveTank1\dataFiles\';
% pathToImageStorage    = '\\Client\H$\ENME489ISpring2012\codes\SPH2DCPPCuda\standardWaveTank1\imageFiles\';
% pathToFffmpeg         = [pathToSPHCode 'ffmpeg01/bin/'];
% pathToHereFromFfmpeg  = ['\\Client\H$\ENME489ISpring2012\codes\DEM2DCPPCuda\standardWaveTank1\imageFiles\'];

pathToSPHCode          = '..\';
pathToDataFileStorage1 = [pathToSPHCode 'dataFiles\'];
pathToDataFileStorage2 = ['X:\Ayan\SPH2DCPPCudaBruteForceGlobalLeapFrog\cases\damBreak_timingvalidation\dataFiles\'];
pathToImageStorage     = [pathToSPHCode 'imageFiles\'];
pathToFffmpeg          = [pathToSPHCode 'ffmpeg01/bin/'];

exportFrames = 0;

A = dir([pathToDataFileStorage1 '*.txt']);
nFiles = length(A);
ind1 = 1;

%get color scaling
plottingAxis = [-0.1497    0.6830   -0.0626    0.5790];
%plottingAxis = [0.4067    0.5592   -0.0057    0.1118];
%plottingAxis = [ 0        1200        1300        2400];
%plottingAxis = [ 0        8800        0        2400];

%for ind1 = [1:1:10]  

while ind1 <= nFiles

   %refresh nFiles
   A = dir([pathToDataFileStorage1 '*.txt']);
   nFiles = length(A);
   s = readInDataFile([pathToDataFileStorage1 A(ind1).name]);
   h = plot02(s,plottingAxis);
   
   %refresh nFiles
   if 1
   A = dir([pathToDataFileStorage2 '*.txt']);
   s = readInDataFile([pathToDataFileStorage2 A(ind1).name]);
   hold on
   h = plot02(s,plottingAxis,'o');
   end
   
   axis equal
   
   
   %pause
      ind1 %loop iteration
      
      
   if exportFrames
   set(gcf,'visible','off','renderer','opengl');
   axis off
   end
   
   %hold on
   %quiver(s(:,1),s(:,2),s(:,3),s(:,4),'b.')
   %set(gca,'dataaspectratio',[1 1 1]);
   %axis([-11   11   -2   15]); - original Funnel axis
   axis([plottingAxis]); % zoom on opening
   
%   colormap copper
%colormap hot
%colormap winter
   C = colormap(jet); %define map
   C = C./max(C(:));
   colormap(C);


   caxis([0 30])
   
   hold off
     
   drawnow
   
   str = A(ind1).name;
   
   if exportFrames;cd(pathToImageStorage);print(gcf,'-dpng',[str(1:end-4) ,'.png']);end
   display([num2str(ind1) ' / ' num2str(nFiles)]);
   ind1 = ind1+1; %this does not affect "for" loop
   %pause
   %keyboard
   if 0
   cd('C:\Users\Chris\Desktop\SPH2DCPPCudaAllPairsv3\postProcessing\')
   save(['v004timeStep'],'s')
   return
   end
   
   
end




%path to here

% %I = regexp(pathToFffmpeg,'/'); %replace / with \
% pathToFffmpeg(I) = '\'; 
% 
% %I = regexp(pathToHereFromFfmpeg,'/'); %replace / with \
% pathToHereFromFfmpeg(I) = '\'; 
% 
%    makeMovieString = [pathToFffmpeg 'ffmpeg -i ' pathToHereFromFfmpeg '\dataFile%05d.png -r 30 ' pathToHereFromFfmpeg '\movieZoom1.mp4'];
%    dos(makeMovieString);
% % % % 
%    makeMovieString2 = [pathToFffmpeg 'ffmpeg -i ' pathToHereFromFfmpeg '\dataFile%05d.png -r 30 -s 800x600 ' pathToHereFromFfmpeg '\movieSmallZoom1.mp4'];
%    dos(makeMovieString2);
% % % 
% % 
% % 
% 




