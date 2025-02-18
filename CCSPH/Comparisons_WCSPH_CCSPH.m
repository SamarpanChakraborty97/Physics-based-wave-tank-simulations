clc
clear all

ele_wcsph = readmatrix("WCSPH.txt");
ele_ccsph = readmatrix("CCSPH.txt");

t_wcsph = readmatrix("Time_WCSPH.txt");
t_ccsph = readmatrix("Time_CCSPH.txt");

ele_wcsph(:,[1,3,5,7,10])=[];
ele_ccsph(:,[1,3,5,7,10])=[];
t_wcsph(:,[1,3,5,7,10])=[];
t_ccsph(:,[1,3,5,7,10])=[];

x_3_experiments = readmatrix("x=3m.csv");
x_3_experiments(end,:) = [];
x_13_experiments = readmatrix("x=13m.csv");
x_22_2_experiments = readmatrix("x=22.2m.csv");
x_32_experiments = readmatrix("x=32m.csv");
x_40_experiments = readmatrix("x=40m.csv");

ele_experiments = horzcat(x_3_experiments(:,2), x_13_experiments(:,2), x_22_2_experiments(:,2), x_32_experiments(:,2), x_40_experiments(:,2));  
t_experiments = horzcat(x_3_experiments(:,1), x_13_experiments(:,1), x_22_2_experiments(:,1), x_32_experiments(:,1), x_40_experiments(:,1));  

Gauges = [3, 13, 22.2, 32, 40];

error_ccsph = zeros(length(Gauges),1);
error_wcsph = zeros(length(Gauges),1);
amp_ratio = zeros(length(Gauges),1);
amp_ratio_c_err = zeros(length(Gauges),1);
amp_ratio_w_err = zeros(length(Gauges),1);

for j=1:length(Gauges)
    loc = Gauges(j);
       
    t = t_experiments(:,j);
    t_start = t(1);
    t_end = t(end);
    index_start = find(t_wcsph(:,j) == t_start);
    index_end = find(t_wcsph(:,j) == t_end);
    
    elev_c = ele_ccsph(index_start:index_end,j);
    elev_w = ele_wcsph(index_start:index_end,j);
    
    mean_ccsph = mean(elev_c);
    mean_wcsph = mean(elev_w);
    mean_experiments = mean(ele_experiments(:,j));

    elev_c = elev_c - mean_ccsph;
    elev_w = elev_w - mean_wcsph;
    elev_exp = ele_experiments(:,j) - mean_experiments;
    
%     ele_5_new = ele_5_read(:,j); 
%     ele_10_new = ele_10_read(:,j);
%     ele_20_new = ele_20_read(:,j);

    figure;
    f = fit(t,elev_c,'smoothingspline','SmoothingParam',0.99);
    h = plot(f);
    xi_c = get(h,'XData');
    yi_c = get(h,'YData');
    yi_c_hilbert = mean(abs(hilbert(yi_c)));
    hold off;
   
    figure;
    f = fit(t,elev_w,'smoothingspline','SmoothingParam',0.99);
    h = plot(f);
    xi_w = get(h,'XData');
    yi_w = get(h,'YData');
    yi_w_hilbert = mean(abs(hilbert(yi_w)));
    hold off;

    figure;
    f2 = fit(t,elev_exp,'smoothingspline','SmoothingParam',0.99);
    h = plot(f2);
    xi_exp = get(h,'XData');
    yi_exp = get(h,'YData');
    yi_exp_hilbert = mean(abs(hilbert(yi_exp)));
    hold off;
 
    error_ccsph(j) = mean(abs(yi_c-yi_exp));
    error_wcsph(j) = mean(abs(yi_w-yi_exp));

    amp_ratio(j) = yi_c_hilbert / yi_w_hilbert;
    amp_ratio_c_err(j) = yi_c_hilbert / yi_exp_hilbert;
    amp_ratio_w_err(j) = yi_w_hilbert / yi_exp_hilbert;
    
    figure;
    g = gcf;
    ax = gca;
    ax.FontSize = 15;
    plot(xi_w, yi_w,'r','Linewidth',2)
    hold on;
    plot(xi_c, yi_c,'b','Linewidth',2)    
    plot(xi_exp, yi_exp,'k','Linewidth',2);
    grid on;
    hold off;
    %xlim([500 800]);
    %ylim([-0.15 0.15]);
    
    name = sprintf("Location = %1.2f m",loc);
    legend('WCSPH','CCSPH','Experimental','FontName','Times','FontSize',12, 'Numcolumns',1, 'Location', 'EastOutside')
    title(name, 'FontName','Times','FontSize',13, 'Fontweight','normal');
    xlabel('Time (seconds)','FontSize', 12)
    ylabel('\eta (metres)','FontSize',12)
    
    g.PaperUnits = 'inches';
    g.PaperPosition = [0 0 16.5 3];
    name2 = sprintf('Comparisons between simulations and experiments at x = %1.2f m_New.png', loc);
    saveas(g,name2)
    close all
end


figure;
hold off;
g = gcf;
ax = gca;
ax.FontSize = 15;
plot(Gauges,amp_ratio,'k-.','Linewidth',2,'marker','*','markersize',6)
grid on;
xlim([0 50]);
xlabel('Gauge locations (in metres)','FontSize', 12)
ylabel('CCSPH wave amplitude/ WCSPH wave amplitude','FontSize',12)

g.PaperUnits = 'inches';
g.PaperPosition = [0 0 7 7];
name2 = sprintf('Comparisons of amplitudes using CCSPH and WCSPH.png');
saveas(g,name2)
close all;

figure;
hold off;
g = gcf;
ax = gca;
ax.FontSize = 15;
plot(Gauges,amp_ratio_c_err,'r-.','Linewidth',2,'marker','s','markersize',6)
hold on;
plot(Gauges,amp_ratio_w_err,'b-.','Linewidth',2,'marker','*','markersize',6)
yline(1,'k--','Linewidth',2);
grid on;
xlim([0 50]);
legend('Using CCSPH','Using WCSPH','FontName','Times','FontSize',12, 'Numcolumns',2, 'Location', 'southwest')
xlabel('Gauge locations (in metres)','FontSize', 12)
ylabel('Experimental wave amplitude/Theoretical wave amplitude','FontSize',12)

g.PaperUnits = 'inches';
g.PaperPosition = [0 0 10 7];
name2 = sprintf('Comparisons of error amplitudes using CCSPH and WCSPH.png');
saveas(g,name2)
close all;

figure;
g = gcf;
ax = gca;
ax.FontSize = 15;
plot(Gauges,error_wcsph,'r-.','Linewidth',2,'marker','s','markersize',6)
hold on;
plot(Gauges,error_ccsph,'b-.','Linewidth',2,'marker','s','markersize',6)
grid on;
xlim([0 45]);
legend('WCSPH','CCSPH','FontName','Times','FontSize',12, 'Numcolumns',1, 'Location', 'northwest')
xlabel('Gauge locations (in metres)','FontSize', 12)
ylabel('Errors','FontSize',12)
g.PaperUnits = 'inches';
g.PaperPosition = [0 0 10 7];
name2 = sprintf('Error comparisons.png');
saveas(g,name2)
