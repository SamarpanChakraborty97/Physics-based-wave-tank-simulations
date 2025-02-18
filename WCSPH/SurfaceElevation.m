close all
clear all

pathToHere         = mfilename('fullpath');
I = regexp(pathToHere,'\');
pathToHere = pathToHere(1:I(end));
pathToDataFileStorage1 = [pathToHere 'dataFiles/'];

storageStride  = 4000;
dt             = 1*10^(-5);

ind = 1;
A = dir([pathToDataFileStorage1 '*.txt']);
nFiles = length(A);
s = readmatrix([pathToDataFileStorage1 A(ind).name]);
I = find(s(:,4)==7);
Freeparticles = s(I,1:2);
dx = 0.01;

Gauges = [0.02 3 6 9 12 15];
L = length(Gauges);

for k=1:L
    I1 = find(Freeparticles(:,1)>Gauges(k)-dx & Freeparticles(:,1)<= Gauges(k)+dx);
    GridFreeParticles = Freeparticles(I1,1:2);
    initial_height = max(GridFreeParticles(:,2))

    ind1=2;
    E1=[];
    T1=[];
    parfor ind1 =1:nFiles
        ind1;
        A = dir([pathToDataFileStorage1 '*.txt']);
        nFiles = length(A);
        s = readmatrix([pathToDataFileStorage1 A(ind1).name]);
        I = find(s(:,4)==7);
        tStep = (ind1-1)*storageStride;
        tTime = tStep*dt;

        Freeparticles = s(I,1:2);
   
        dx = 0.01;
        I1 = find(Freeparticles(:,1)>Gauges(k)-dx & Freeparticles(:,1)<= Gauges(k)+dx);    
        GridFreeParticles = Freeparticles(I1,1:2);
        if ~isempty(GridFreeParticles(:,1))
            surface_elevation = max(GridFreeParticles(:,2)) - initial_height;
        else
            surface_elevation = 0;
        end
        E1(ind1) = surface_elevation;
        T1(ind1) = tTime;
   
    %ind1 = ind1+1;
    end
    Elevation1 = E1';
    Time1 = T1';

    f1 = fit(Time1,Elevation1,'smoothingspline','SmoothingParam',0.999);
    h1 = plot(f1);
    xi1 = get(h1,'XData');
    yi1 = get(h1,'YData');

    data=[xi1' yi1'];
    filename = sprintf('Case1_x=%d.txt',Gauges(k));
    dlmwrite(filename,data);
end