close all
clear all

pathToHere         = mfilename('fullpath');
I = regexp(pathToHere,'\');
pathToHere = pathToHere(1:I(end));
pathToDataFileStorage = [pathToHere 'dataFiles/'];

storageStride  = 4000;
dt             = 1*10^(-5);
% height  = 0.3;
A  = dir([pathToDataFileStorage '*.txt']);
nFiles = length(A);
S = 0.1657;
H = 0.15;
a = 0.0747;
T = 2.0;
L = 4.52;
k = 2*pi/L;
w = 2*pi/T;
d = 0.66;

ind = 1;
A = dir([pathToDataFileStorage '*.txt']);
nFiles = length(A);
s = readInDataFile([pathToDataFileStorage A(ind).name]);
I = find(s(:,7)==7);
Freeparticles = s(I,1:2);
dx = 0.1;
I1 = find(Freeparticles(:,1)>3.5 -dx & Freeparticles(:,1)<= 3.5+dx);
GridFreeParticles = Freeparticles(I1,1:2);
initial_height = max(GridFreeParticles(:,2))

ind1=2;
E=[];
T=[];
while ind1 <= nFiles
    ind1;
    A = dir([pathToDataFileStorage '*.txt']);
    nFiles = length(A);
    
   s = readInDataFile([pathToDataFileStorage A(ind1).name]);
   
   tStep = (ind1-1)*storageStride;
   tTime = tStep*dt;
   
   I = find(s(:,7)==7);
   Freeparticles = s(I,1:2);
   
   dx = 0.1;
   I1 = find(Freeparticles(:,1)>3.5 -dx & Freeparticles(:,1)<= 3.5+dx);    
   GridFreeParticles = Freeparticles(I1,1:2);
   surface_elevation = max(GridFreeParticles(:,2)) - initial_height;
   E(ind1) = surface_elevation
   T(ind1) = tTime;
   
   ind1 = ind1+1;
end
Elevation = E';
Time = T';
f = fit(Time,Elevation,'smoothingspline','SmoothingParam',0.99);
h = plot(f);
xi = get(h,'XData');
yi = get(h,'YData');

t = linspace(0,10,500);
% x  = 6.5;
%ele = a * sin(k*x - w*t) - (k*a^2 /2)*((3 - (tanh(k * d))^2)/(2 * (tanh(k * d))^2)) * cos(2*k*x - 2 *w*t) - k*a^2/(2 *sinh(2*k*d));
% ele = a * cos(w*t - k*x);
% plot(t , ele ,'r');
% hold on
plot(xi,yi,'k') ;
grid on;
% legend('Theoretical Stokes 1st order','Simulation results')
xlabel('Time (seconds)')
ylabel('\eta (metres)')
% title('Surface elevation comparison between theory and simulation results')
f = gcf;
saveas(f,'Elevation_latest2.png')
   