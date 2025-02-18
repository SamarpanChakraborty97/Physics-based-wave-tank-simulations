clear all
clc

a = readmatrix("spectrum076.txt");

Tp = 12.50;
Hs = 2.20;
d = 27.4;

Nf = length(a(:,1));
Freqs = a(:,1);
Energy = a(:,3);

del_f = a(1,2);

amps = sqrt(2.*Energy*del_f);
Omegas = 2*pi.*Freqs;

K = zeros(Nf,1);
for i=1:Nf
    syms ki
    eqn = (2*pi*Freqs(i))^2 == 9.81*ki*tanh(ki*27.4);
    k = vpasolve(eqn);
    K(i) = k;
end

k_p = zeros(1,1);
fp = 1/Tp;
wp = 2*pi*fp;
syms kp;
assume (kp > 0);
eqn = (wp^2 == 9.81*kp*tanh(kp*d));
k_p(1) = vpasolve(eqn);
lambda_p = abs(round(2*pi/k_p, 2));

t = 0;
dt = 1 * 10^(-3);
T = [];
P = [];
x = -50;
while t<200
    sum = 0;
    for i = 1:Nf
        k = K(i);
        amp = amps(i);
        w = Omegas(i);
        sum = sum + amp * cos(k*x - w*t);
    end
    T(end+1) = t;
    P(end+1) = sum;
    t = t + dt;
end

plot(T,P,'r','Linewidth',2)
xlim([T(1) T(end)]);
grid on;
grid minor;
xlabel('$$t(secs)$$','interpreter','latex')
ylabel('$$S(t) (metres)$$','interpreter','latex')
title('Piston movement time series','interpreter','latex')
f = gcf;
f.PaperUnits = 'inches';
f.PaperPosition = [0 0 10 3];
set(gca,'FontName','Times','FontSize',12);
print(f,'-dpng','Piston movement time series.png')
saveas(f,'Piston movement time series2.png')
