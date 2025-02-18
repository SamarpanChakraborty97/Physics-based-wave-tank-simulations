clc
clear all

ele_ccsph = readmatrix("Elevation_CCSPH_shallowWater.txt");
ele_non_ccsph = readmatrix("Elevation_Leapfrog_shallowWater.txt");

time = readmatrix("Time_CCSPH_shallowWater.txt");

S = 0.1657;
H = 0.15;
a = 0.075;
T = 2.0;
L = 4.52;
k = 2*pi/L;
w = 2*pi/T;
d = 0.66;

num_gauges = length(ele_ccsph(1,:));
num_entries = length(ele_ccsph(:,1));

Gauges = [0.5, 2, 4, 6, 8, 10, 12, 14, 15, 20, 25, 30, 35, 40, 45];

amp_ratio = zeros(length(Gauges),1);
amp_c_ratio = zeros(length(Gauges),1);
amp_n_ratio = zeros(length(Gauges),1);

error_c = zeros(length(Gauges),1);
error_n = zeros(length(Gauges),1);

for j=1:num_gauges
       
    t = time(:,j);
    mean_c = mean(ele_ccsph(:,j));
    mean_n = mean(ele_non_ccsph(:,j));
    
    ele_c = ele_ccsph(:,j)-mean_c;
    ele_n = ele_non_ccsph(:,j)-mean_n;
    
    %ele_c = ele_ccsph(:,j);
    %ele_n = ele_non_ccsph(:,j);
    
    loc = Gauges(j);
    ele_theo = a*(cos(k*loc - w*t)+k*a*((3-tanh(k*d)*tanh(k*d))/(4*tanh(k*d)*tanh(k*d)*tanh(k*d)))*cos(2*(k*loc - w*t)));
    %ele_theo = a * sin(k*loc - w*t) - (k*a^2 /2)*((3 - (tanh(k * d))^2)/(2 * (tanh(k * d))^2)) * cos(2*k*loc - 2 *w*t) - k*a^2/(2 *sinh(2*k*d));

    f = fit(t,ele_c,'smoothingspline','SmoothingParam',0.99);
    h = plot(f);
    xi_c = get(h,'XData');
    yi_c = get(h,'YData');
    
    index = find(xi_c >= 74.8 & xi_c < 75.2);
    ratio = floor(index(1)/75);   
    
    yi_c_req_hilbert = mean(abs(hilbert(yi_c(floor(ratio*(loc+1)):index))));
      
    f2 = fit(t,ele_n,'smoothingspline','SmoothingParam',0.99);
    h2 = plot(f2);
    xi_n = get(h2,'XData');
    yi_n = get(h2,'YData');
    
    index = find(xi_n >= 74.8 & xi_n < 75.2);
    ratio = floor(index(1)/75);   
    yi_n_req_hilbert = mean(abs(hilbert(yi_n(floor(ratio*(loc+1)):index))));
    
    f3 = fit(t,ele_theo,'smoothingspline','SmoothingParam',0.99);
    h3 = plot(f3);
    xi_t = get(h3,'XData');
    yi_t = get(h3,'YData');
    
    index = find(xi_t >= 74.8 & xi_t < 75.2);
    ratio = floor(index(1)/75);   
    yi_t_req_hilbert = mean(abs(hilbert(yi_t(floor(ratio*(loc+1)):index))));
    
    error_c(j) = mean(abs(ele_c-ele_theo));
    error_n(j) = mean(abs(ele_n-ele_theo));
    
    amp_ratio(j) = yi_c_req_hilbert / yi_n_req_hilbert;
    amp_c_ratio(j) = yi_c_req_hilbert / yi_t_req_hilbert;
    amp_n_ratio(j) = yi_n_req_hilbert / yi_t_req_hilbert;
   
    figure;
    g = gcf;
    ax = gca;
    ax.FontSize = 15;
    plot(xi_c,yi_c,'k','Linewidth',2)
    hold on;
    plot(xi_n,yi_n,'r','Linewidth',2);
    plot(t,ele_theo,'b','Linewidth',2);
    grid on;
    hold off;
    xlim([40 75]);
    
    
    name = sprintf("Location = %1.2f m",loc);
    legend('Using CCSPH','Without using CCSPH','Theoretical results','FontName','Times','FontSize',12, 'Numcolumns',1, 'Location', 'EastOutside')
    title(name, 'FontName','Times','FontSize',13, 'Fontweight','normal');
    xlabel('Time (seconds)','FontSize', 12)
    ylabel('\eta (metres)','FontSize',12)
    
    g.PaperUnits = 'inches';
    g.PaperPosition = [0 0 16.5 3];
    name2 = sprintf('Comparisons between CCSPH and normal at x = %1.2f m for shallow waters_LongerTank_New.png', loc);
    saveas(g,name2)
    
    hold off;
    
end

figure;
hold off;
g = gcf;
ax = gca;
ax.FontSize = 15;
plot(Gauges,amp_ratio,'k-.','Linewidth',2,'marker','*','markersize',6)
grid on;
xlim([0 50]);
%legend('Amplitude ratios of method using CCSPH and WCSPH over wave tank','FontName','Times','FontSize',9, 'Numcolumns',2, 'Location', 'northwest')
xlabel('Gauge locations (in metres)','FontSize', 12)
ylabel('CCSPH wave amplitude / WCSPH wave amplitude','FontSize',12)

g.PaperUnits = 'inches';
g.PaperPosition = [0 0 7 7];
name2 = sprintf('Comparisons of amplitudes using CCSPH and not for shallow waters_LongerTank_New.png');
saveas(g,name2)
close all;

figure;
hold off;
g = gcf;
ax = gca;
ax.FontSize = 15;
plot(Gauges,amp_c_ratio,'r-.','Linewidth',2,'marker','s','markersize',6)
hold on;
plot(Gauges,amp_n_ratio,'b-.','Linewidth',2,'marker','*','markersize',6)
yline(1.00, 'k-.','Linewidth',2);
grid on;
xlim([0 50]);
legend('Using CCSPH','Not using CCSPH','FontName','Times','FontSize',12, 'Numcolumns',2, 'Location', 'southwest')
xlabel('Gauge locations (in metres)','FontSize', 12)
ylabel('Experimental wave amplitude / Theoretical wave amplitude','FontSize',9)

g.PaperUnits = 'inches';
g.PaperPosition = [0 0 10 7];
name2 = sprintf('Comparisons of amplitude ratios for shallow waters_LongerTank_New.png');
saveas(g,name2)
close all;

figure;
hold off;
g = gcf;
ax = gca;
ax.FontSize = 15;
plot(Gauges,error_c,'r-.','Linewidth',2,'marker','s','markersize',6)
hold on;
plot(Gauges,error_n,'b-.','Linewidth',2,'marker','*','markersize',6)
grid on;
xlim([0 50]);
legend('Using CCSPH','Not using CCSPH','FontName','Times','FontSize',12, 'Numcolumns',2, 'Location', 'northwest')
xlabel('Gauge locations (in metres)','FontSize', 12)
ylabel('errors (in metres)','FontSize',12)

g.PaperUnits = 'inches';
g.PaperPosition = [0 0 10 7];
name2 = sprintf('Comparisons of errors for shallow waters_LongerTank_New.png');
saveas(g,name2)
close all;
    
