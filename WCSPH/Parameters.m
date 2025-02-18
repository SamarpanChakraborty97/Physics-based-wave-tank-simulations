clc
clear all

fl = 0.8;
fh = 4.0;
fp = 1.82;
Tp = 1/fp;
H_13 = 0.0075;
g = 4;
h = 0.2;
alpha = (0.0624/(0.230 + 0.0336*g - 0.185/(1.9+g)));
Nf = 80;
delf = (fh - fl)/Nf;
Ab = 0.00375;
F = zeros(Nf,1);
Sf = zeros(Nf,1);
% tb = 26.8328;
xb = 3.5;

F(1) = fl;
Sf(1) = alpha * (H_13^2)*(Tp^(-4)*(F(1))^(-5))*exp(-1.25*(Tp*F(1))^-4)*(g^(exp(-(F(1)/fp)-1)^2 / (2*0.07^2)));

for i=2:Nf
    F(i) = F(i-1) + delf;
    if F(i)>=fp
        s = 0.09;
        Sf(i) = alpha * (H_13^2)*(Tp^(-4)*(F(i))^(-5))*exp(-1.25*(Tp*F(i))^-4)*(g^(exp(-(F(i)/fp)-1)^2 / (2*s^2)));
    else
        s = 0.07;
        Sf(i) = alpha * (H_13^2)*(Tp^(-4)*(F(i))^(-5))*exp(-1.25*(Tp*F(i))^-4)*(g^(exp(-(F(i)/fp)-1)^2 / (2*s^2)));
    end
end

% a =  Ab * Sf/ sum(Sf);
%S = sum(Sf * delf);
for i=1:Nf
    syms ki
    eqn = (2*pi*F(i))^2 == 9.81*ki*tanh(ki*0.2);
    k = vpasolve(eqn);
    K(i) = k;
end

a =  sqrt(2*Sf*delf);

for i=1:Nf
    B(i) = 2*(sinh(K(i)*h)*sinh(K(i)*h))/(sinh(K(i)*h)*cosh(K(i)*h)+ K(i)*h);
% B = 2*(sinh(K*h)*sinh(K*h))/(sinh(K*h)*cosh(K*h))+K*h
end

P(:,1) = a;
P(:,2) = B;
P(:,3) = F;
P(:,4) = K';
P(:,5) = Sf;

writematrix(P,'InitialParameters.txt')

% t=0:1e-5:40;
% Eta = zeros(length(t),1);
% for i = 1:length(t)
%     sum = 0;
%     time = t(i);
%     for j=1:Nf
%         sum = sum + a(j)*cos(K(j)*(-xb) - 2*pi*F(j)*(time-tb));
%     end
%     Eta(i) = sum;
% end
% plot(t,Eta)
% 

 
% for i=1:Nf
%     syms ki
%     eqn = (2*pi*F(i))^2 == 9.81*ki*tanh(ki*1.2);
%     vpasolve(eqn,ki)
%     K(i) = k;
% end
% K