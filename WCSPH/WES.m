clc
clear all

Gauges = [2.000000e-02,3,6,9,12,15];

name = sprintf('Case7_x=%d.txt',Gauges(3));
m = readmatrix(name);
Ele = m(:,2);
T = m(;,1);

Fb = 2;
Fc = 1.75;

freq = 0:3/70:3;
a = freq * Fc;
lF = length(a);

dTau = 4/7;
Tau = 0:dTau:40;
lT = length(Tau);

WES = zeros(lF,lT);
for i = 1:length(Tau)
    for j = 1:length(a)
        sum = 0;
        for k = 1:length(Ele)
            sum = sum + Ele(k) * conj((1/sqrt(pi * Fb)) * exp(2*pi *1i * Fc * (T(k) - Tau(i))/a(j)) * exp(-((T(k) - Tau(i))/a(j))^2)/Fb);
        end
        Xw = 1/sqrt(a(j)) * sum;
        WES(i,j) = Xw * conj(Xw) / a(j);
    end
end
WES



    
    
