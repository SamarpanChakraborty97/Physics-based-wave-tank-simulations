clc
clear all

matrix = readmatrix('Phases.txt');
phaseValues1 = matrix(:,1);
phaseValues2 = matrix(:,2);
phaseValues3 = matrix(:,3);
phaseValues4 = matrix(:,4);
phaseValues5 = matrix(:,5);
phaseValues6 = matrix(:,6);
% length(phaseValues1(:,1))

% plot(phaseValues1)
% hold on;
% plot(phaseValues2)
% plot(phaseValues3)
% plot(phaseValues4)
% plot(phaseValues5)
% plot(phaseValues6)

a0 = 0.0281;
k0 = 0.8;
w0 = 2.97;
a = 0.784;
eI = k0 * a0;
s2 = sqrt(2);
storageStride = 4000;
dt = 0.00001;


x = 0;
X = [];
X(1) = x;
T = [];
t=0;
T(1) = t;
% (eI*eI*eI/(4*s2))*w0*t


% sin(phaseValues1(counter))
% a0 * sech((eI*s2)*phaseValues1(counter) - (eI*eI*eI/(4*s2))*w0*t) * cos(phaseValues1(counter));
while t<40
    counter = floor(t/(storageStride*100*dt))+1;
    w1 = a0 * 1/cosh((eI*s2)*phaseValues1(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(phaseValues1(counter));
    w2 = a0 * 1/cosh((eI*s2)*phaseValues2(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(phaseValues2(counter));
    w3 = a0 * 1/cosh((eI*s2)*phaseValues3(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(phaseValues3(counter));
    w4 = a0 * 1/cosh((eI*s2)*phaseValues4(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(phaseValues4(counter));
    w5 = a0 * 1/cosh((eI*s2)*phaseValues5(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(phaseValues5(counter));
    w6 = a0 * 1/cosh((eI*s2)*phaseValues6(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(phaseValues6(counter));
    
    x = a * (w1+w2+w3+w4+w5+w6);
    t = t + dt;
    T(end+1) = t;
    X(end+1) = x;
end
figure
plot(T,X)
% xlim([0 5e-2])
