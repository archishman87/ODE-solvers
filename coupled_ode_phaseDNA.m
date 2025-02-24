function [] = coupled_ode_phaseDNA(d, DNAmax)
    % Constants
    T7 = 470;
    
    
    % Parameters
    kr = 0.016079516;
    kb1 = 1.22851158;
    kb2 = 0.000100154;
    kp = 0.107670335;
    k1 = 22.91161736;
    k2 = 1.958401957;
    k3 = 1.956331343;
    kmat = 1.957970722;
    kd = 0.004291751;
  
    TsR = 91.69390911;
    TlR = 905.0532033;
    
    
    params = [kr, kb1, kb2, kp, k1, k2, k3, kmat, kd];
    
    tspan = [0 6];
    
    % Generate a colormap with jet scheme
    %figure;
    hold all;
    DNAall = [0:d:DNAmax];
    colorMap = jet(length(DNAall));
    
    % Plot trajectories and calculate arrows
    for j = 1:length(DNAall)
        DNA = DNAall(j);
        RNase = 370;

        initial_conditions = [0, 0, 0, 0, 0, 0, 0, TsR, TlR];
        [t, y] = ode45(@(t,y) coupled_odes(t, y, params, [T7, DNA, RNase]), tspan, initial_conditions);

        % Plot trajectory

        color = colorMap(j, :);
        %plot(t,y(:,2),'Color', color,'LineWidth',2,'LineStyle','-');
        plot(t,y(:,7),'Color', color,'LineWidth',2,'LineStyle','-');
        % g = y(:,7);
        % m = y(:,2);
        % [maxYValue, indexAtMaxY] = max(y(:,2));
        % xValueAtMaxYValue = g(indexAtMaxY(1));
        % plot(xValueAtMaxYValue,maxYValue,'k.','MarkerSize',10)

        % Calculate arrows based on trajectory points
        % dt = 10; % Time step for plotting arrows (adjust as needed)
        % for k = 1:dt:length(t)-1
        %     % Calculate direction at each point along the trajectory
        %     dydt = coupled_odes(t(k), y(k,:), params, constants);
        %     u = dydt(7);
        %     v = dydt(2);
        % 
        %     % Plot arrow at current point
        %     quiver(y(k,7), y(k,2), u, v, 'Color', color, 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
        % end
        
    end
    
    % h = colorbar;
    % %set(h, 'ylim', [1 length(TsRall)])
    % ticks = linspace(0, 1, 5);
    % ticklabels = {'0','500','1000','1500','2000'};
    % h.Ticks = ticks;
    % h.TickLabels = ticklabels;
    % 
    % axis tight;
    % colormap(jet);  % Set the colormap to match our lines
    % xlabel('Time (hours)');
    % ylabel('[mRNA]/[GdmS^*] (nM)');
    % %title('Phase Plot with Time Arrows');
    % hold off;
    % grid on
    % set(gca,'FontSize',15,'xcolor','k','ycolor','k','FontWeight','bold')
    % ax = gca;
    % ax.YColor = 'g';
    % ax.XColor = 'm';
end
