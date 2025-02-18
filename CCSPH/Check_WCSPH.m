ele_wcsph = readmatrix("WCSPH.txt");
t_wcsph = readmatrix("Time_WCSPH.txt");

ele_wcsph(:,[1,3,5,7,10])=[];
t_wcsph(:,[1,3,5,7,10])=[];

Gauges = [3, 13, 22.2, 32, 40];

for j=1:length(Gauges)
    loc = Gauges(j);
    
    elev_w = ele_wcsph(:,j);
    mean_wcsph = mean(elev_w);
    
    t_w = t_wcsph(:,j);
    
    figure;
    f = fit(t_w,elev_w,'smoothingspline','SmoothingParam',0.99);
    h = plot(f);
    xi_w = get(h,'XData');
    yi_w = get(h,'YData');
    yi_w_hilbert = mean(abs(hilbert(yi_w)));
    hold off;
    
    figure;
    g = gcf;
    ax = gca;
    ax.FontSize = 15;
    plot(xi_w, yi_w*100,'r','Linewidth',2)
    name = sprintf("Location = %1.2f m",loc);
    %legend('WCSPH','CCSPH','Experimental','FontName','Times','FontSize',12, 'Numcolumns',1, 'Location', 'EastOutside')
    title(name, 'FontName','Times','FontSize',13, 'Fontweight','normal');
    xlabel('Time (seconds)','FontSize', 12)
    ylabel('\eta (metres)','FontSize',12)
    xlim([20 50])
    ylim([-6 6])
    
    g.PaperUnits = 'inches';
    g.PaperPosition = [0 0 16.5 3];
    name2 = sprintf('Check for WCSPH x = %1.2f m_New.png', loc);
    saveas(g,name2)
    close all
end
