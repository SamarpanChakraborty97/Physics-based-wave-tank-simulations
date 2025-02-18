P = readmatrix('InitialParameters.txt');

a = P(:,1);
F = P(:,3);
K = P(:,4);
Nf = 80;
xb = 3.0;
tb = 15;

t=0:0.01:40;
Eta = zeros(length(t),1);
for i = 1:length(t)
    sum = 0;
    time = t(i);
    for j=1:Nf
        sum = sum + a(j)*cos(K(j)*(xb) - 2*pi*F(j)*(time-tb));
    end
    Eta(i) = sum;
end
plot(t,Eta*100,'k','Linewidth',2)
% xlim([0 60]);
grid on;
grid minor;
f = gcf;
ax = gca;
set(ax,'FontName','Times')
f.Position = [100 100 1200 300];
ylabel('X (cm)');
xlabel('Time (seconds)');
