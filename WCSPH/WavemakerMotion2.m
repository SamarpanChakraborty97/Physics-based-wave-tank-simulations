clc
clear all

matrix = readmatrix('Phases.txt');
phaseValues1 = matrix(:,1);
phaseValues2 = matrix(:,2);
phaseValues3 = matrix(:,3);
phaseValues4 = matrix(:,4);
phaseValues5 = matrix(:,5);
phaseValues6 = matrix(:,6);

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

alpha = a0 * k0 * w0/2;
beta = a0^2*k0^2*w0/4 + w0;

while t<40
    counter = floor(t/(storageStride* 100 *dt))+1;
    w1 = a0 * 1/cosh(phaseValues1(counter) - alpha*t)*cos(beta * t + phaseValues1(counter));
    w2 = a0 * 1/cosh(phaseValues2(counter) - alpha*t)*cos(beta * t + phaseValues2(counter));
    w3 = a0 * 1/cosh(phaseValues3(counter) - alpha*t)*cos(beta * t + phaseValues3(counter));
    w4 = a0 * 1/cosh(phaseValues4(counter) - alpha*t)*cos(beta * t + phaseValues4(counter));
    w5 = a0 * 1/cosh(phaseValues5(counter) - alpha*t)*cos(beta * t + phaseValues5(counter));
    w6 = a0 * 1/cosh(phaseValues6(counter) - alpha*t)*cos(beta * t + phaseValues6(counter));
    
    x = a * (w1+w2+w3+w4+w5+w6);
    t = t + dt;
    T(end+1) = t;
    X(end+1) = x;
end
figure
plot(T,X)