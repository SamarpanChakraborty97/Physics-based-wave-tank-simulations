clc
clear all

matrix = readmatrix('Phases.txt');
phaseValues1 = matrix(1,:);
phaseValues2 = matrix(2,:);
phaseValues3 = matrix(3,:);
phaseValues4 = matrix(4,:);
phaseValues5 = matrix(5,:);
phaseValues6 = matrix(6,:);
% phaseValues7 = matrix(7,:);
% phaseValues8 = matrix(8,:);
% phaseValues9 = matrix(9,:);
% phaseValues10 = matrix(10,:);
% phaseValues11 = matrix(11,:);
% phaseValues12 = matrix(12,:);
% phaseValues13 = matrix(13,:);
% phaseValues14 = matrix(14,:);
% phaseValues15 = matrix(15,:);
% length(phaseValues1(:,1))

% plot(phaseValues1)
% hold on;
% plot(phaseValues2)
% plot(phaseValues3)
% plot(phaseValues4)
% plot(phaseValues5)
% plot(phaseValues6)

a0 = 0.055;
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
beta = a0^2*k0^2*w0/4 + w0;

% sin(phaseValues1(counter))
% a0 * sech((eI*s2)*phaseValues1(counter) - (eI*eI*eI/(4*s2))*w0*t) * cos(phaseValues1(counter));
while t<40
    counter = floor(t/(storageStride*100*dt))+1;
    w1 = a0 * 1/cosh((eI*s2)*phaseValues1(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues1(counter));
    w2 = a0 * 1/cosh((eI*s2)*phaseValues2(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues2(counter));
    w3 = a0 * 1/cosh((eI*s2)*phaseValues3(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues3(counter));
    w4 = a0 * 1/cosh((eI*s2)*phaseValues4(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues4(counter));
    w5 = a0 * 1/cosh((eI*s2)*phaseValues5(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues5(counter));
    w6 = a0 * 1/cosh((eI*s2)*phaseValues6(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues6(counter));
%     w7 = a0 * 1/cosh((eI*s2)*phaseValues7(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues7(counter));
%     w8 = a0 * 1/cosh((eI*s2)*phaseValues8(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues8(counter));
%     w9 = a0 * 1/cosh((eI*s2)*phaseValues9(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues9(counter));
%     w10 = a0 * 1/cosh((eI*s2)*phaseValues10(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues10(counter));
%     w11 = a0 * 1/cosh((eI*s2)*phaseValues11(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues11(counter));
%     w12 = a0 * 1/cosh((eI*s2)*phaseValues12(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues12(counter));
%     w13 = a0 * 1/cosh((eI*s2)*phaseValues13(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues13(counter));
%     w14 = a0 * 1/cosh((eI*s2)*phaseValues14(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues14(counter));
%     w15 = a0 * 1/cosh((eI*s2)*phaseValues15(counter) - (eI*eI*eI/(4*s2))*w0*t)*cos(beta * t + phaseValues15(counter));
    
    x = a * (w1+w2+w3+w4+w5+w6);
    t = t + dt;
    T(end+1) = t;
    X(end+1) = x;
end
f = figure;
plot(T,X);
grid on;
grid minor;
print(f,'PhaseSync1.png','-dpng')

% xlim([0 5e-2])
