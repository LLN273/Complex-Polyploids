function [  totalSumsofDistsToClusterCentroids_1, delta_totalSumsofDistsToClusterCentroids  ] = clustering_control( My_Data, N_runs, Kmax, plot_id, title_info )

% manages clustering analysis & prepares associated plots
% My_Data:          input data array 
% N_runs:           number of runs for each k-value
% Kmax:             maximum k-value evaluated
% plot_id:          plot number



[ totalSumsofDistsToClusterCentroids_1 ] = clustering_analysis_corr_v1( My_Data, N_runs, Kmax );

totalSumsofDistsToClusterCentroids_forward_1 = totalSumsofDistsToClusterCentroids_1; 
totalSumsofDistsToClusterCentroids_forward_2 = totalSumsofDistsToClusterCentroids_1;
totalSumsofDistsToClusterCentroids_forward_1(10,:) =[];
totalSumsofDistsToClusterCentroids_forward_2(1,:) = [];
delta_totalSumsofDistsToClusterCentroids = totalSumsofDistsToClusterCentroids_forward_1 - totalSumsofDistsToClusterCentroids_forward_2;  % drop in totalSumsofDistsToClusterCentroids


% plots
hFig1 = figure(plot_id);
scrsz = get(groot,'ScreenSize');
set(hFig1, 'Position', [1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
subplot(2,1,1)
boxplot(totalSumsofDistsToClusterCentroids_1', 'Labels',{'1', '2','3','4','5','6','7','8','9','10'})
xlabel('k','FontSize',16);
ylabel('S_D','FontSize',16);
set(gca,'FontSize',16)
title(title_info);
subplot(2,1,2)
boxplot(delta_totalSumsofDistsToClusterCentroids', 'Labels',{'k1 to k2','k2 to k3','k3 to k4','k4 to k5','k5 to k6','k6 to k7','k7 to k8','k8 to k9','k9 to k10'})
xlabel('DELTA k','FontSize',16);
ylabel('DELTA S_D','FontSize',16);
set(gca,'FontSize',16)

n_pause = 2;
pause('on')
pause(n_pause)                                               

end

