 clc
    clear all

    h = 2.0;
    Tp = 2.0;
    Hs = 0.6;
    gamma = 3.3;
    alpha = 0.0624 / (0.23 + 0.0336 * gamma - (0.185/(1.9+gamma)));

    fp = (1/Tp);
    fStart = 0.2;
    fStop = 2.0;

    N = 128;
    delF = (fStop - fStart)/ N;
    f = zeros(N,1);
    f2 = f;
    k = f; L = f; S = f; w = f;B = f;H = f;
    for i = 1:N
        f(i) = fStart + i*delF - delF/2;
        f2(i) = f(i)/fp;
        w(i) = 2* pi * f(i);
        k(i) = w(i) * w(i)/9.81;
        L(i) = 2*pi/ k(i);
        B(i) = 2 * sinh(k(i)*h) * sinh(k(i)*h)/((sinh(k(i)*h) * cosh(k(i)*h))+k(i)*h);
        if f(i) < fp
            sigma = 0.07;
            beta = exp(-(f(i)-fp)^2/(2*sigma^2*fp^2));
            S(i) = alpha * Hs^2 * fp^4 * f(i)^(-5) * gamma^beta * exp((-5/4)*(fp/f(i))^4);
        else
            sigma = 0.09;
            beta = exp(-(f(i)-fp)^2/(2*sigma^2*fp^2));
            S(i) = alpha * Hs^2 * fp^4 * f(i)^(-5) * gamma^beta * exp((-5/4)*(fp/f(i))^4);
        end
        H(i) = 2 * sqrt(2 * S(i) * delF);
    end
    del = (2*pi).*rand(N,1);
    t = 0;
dt = 1 * 10^(-5);
T = [];
P = [];
while t<30
    sum = 0;
    for i = 1:N
        sum = sum + H(i)/(2*B(i)) * cos(w(i)*t + del(i)+pi/2);
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
saveas(f,'Piston movement Ssytem 1.png')