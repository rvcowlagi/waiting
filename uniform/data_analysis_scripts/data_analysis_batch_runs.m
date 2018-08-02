% Summarize various runs with randomized paramters including
% 1) Random fields (# of peaks, locations, spreads, weights)
% 2) Rancom Costs  (field exposure scaling, waiting, movement)
clear all; close all; clc;

% load wait_study02_CPORT_fixedField.mat
% load wait_study02_HL_allRandWeights.mat
% load wait_study02_CPORT_250_randFC.mat
% load wait_study02_500randFC.mat
% load wait_study02_randField_zeroWait.mat
% load wait_study02_3peaks_zeroWait.mat
% load wait_study02_VaryingLocations.mat
% load wait_study02_VaryingWidths.mat
% load wait_study02_VaryLocationWidthFmin.mat
% load wait_study02_VaryingLocWidthFmin_N200.mat
% load wait_study02_VaryLocWidthFmin_N500.mat
% load wait_study02_VaryLocWidthFmin_ZeroWait.mat
% load wait_study02_VaryLocWidthFmin_N2000.mat
% load wait_study02_VaryLocWidthFmin_ZeroWait_N1000.mat
% load wait_study02_VaryLocWidthWait_ZeroFmin_N1000.mat
% load wait_study02_VaryLocWidthWait_ZeroFmin_30Peaks.mat
% load wait_study02_VaryLocWidthWaitFmin_40peaks.mat
% load wait_study02_VaryLocWidthWaitFmin_40pks_CompTime.mat
load wait_study02_VaryLocWidthWaitFmin_40pk_LcheckFix.mat

if (exist('common_search_data','var'))
    search_data = common_search_data;
end

n_sim = size(sim_data_record,2);
t_final = search_data.t_max;
n_grid_pt = sqrt(search_data.n_vertices);
wksp = 10;
grid_sep = search_data.grid_sep;
x_grid_srch	= linspace(-wksp + grid_sep/2, wksp - grid_sep/2, n_grid_pt);
y_grid_srch	= linspace(-wksp + grid_sep/2, wksp - grid_sep/2, n_grid_pt);
[x_plot_grid, y_plot_grid] = meshgrid(x_grid_srch, y_grid_srch);

cost_diff = zeros(n_sim,1); 
wait_weight = zeros(n_sim,1); move_weight = zeros(n_sim,1);
exposure_weight = zeros(n_sim,1);
min_difft = zeros(n_sim,1); max_difft = zeros(n_sim,1); min_field = zeros(n_sim,1);

for m1 = 1:n_sim
    fprintf('%i/%i\n',m1,n_sim)
    cost_diff(m1) = sim_data_record(m1).costs.nowait - sim_data_record(m1).costs.wait;
    cost_diff_pcent(m1) = cost_diff(m1)/sim_data_record(m1).costs.nowait*100;
    cost_diff_check(m1) = sim_data_record(m1).costs.nowait - sim_data_record(m1).costs.wait_check;
    cost_diff_check_pcent(m1) = cost_diff_check(m1)/sim_data_record(m1).costs.nowait*100;
    wait_weight(m1) = sim_data_record(m1).weights.wait;
    move_weight(m1) = sim_data_record(m1).weights.move;
    exposure_weight(m1) = sim_data_record(m1).weights.exposure;
    
    % Obtain field parameters
    num_peaks(m1) = sim_data_record(m1).coeffs.n_peaks;
    n_peaks(m1) = sim_data_record(m1).coeffs.n_peaks;
    coeff_peaks_0 = sim_data_record(m1).coeffs.coeff_peaks_0;
    coeff_peaks_f = sim_data_record(m1).coeffs.coeff_peaks_f;
    coeff_peaks_rate= (coeff_peaks_f - coeff_peaks_0)/t_final;
    
    max_Wndot(m1) = max(coeff_peaks_rate(1,:));
    max_xndot(m1) = max(coeff_peaks_rate(2,:));
    max_yndot(m1) = max(coeff_peaks_rate(3,:));
    max_sigxndot(m1) = max(coeff_peaks_rate(4,:));
    max_sigyndot(m1) = max(coeff_peaks_rate(5,:));
    
    field_params.min_assigned = sim_data_record(m1).coeffs.min_assigned;
    field_params.max_assigned = sim_data_record(m1).coeffs.max_assigned;
    field_params.t_final = t_final;
    field_params.n_peaks = n_peaks(m1);
    field_params.coeff_peaks_0 = coeff_peaks_0;
    field_params.coeff_peaks_rate = coeff_peaks_rate;
	tvspf	= threat_field_args(x_plot_grid, y_plot_grid, search_data.t_grid, 'peaks', 0, field_params);
    max_difft(m1) = max(tvspf.difft(:));
    min_difft(m1) = min(tvspf.difft(:));
    min_field(m1) = min(tvspf.field(:));
    max_field(m1) = max(tvspf.field(:));
    avg_field(m1) = mean(tvspf.field(:));
    avg_difft(m1) = mean(tvspf.difft(:));
    clear tvspf;
    
    % Path properties
    path_length_wait(m1)   = length(sim_data_record(m1).paths.wait.v);
    path_length_nowait(m1) = length(sim_data_record(1).paths.nowait);
    
    exec_data_wait_n_iter(m1) = sim_data_record(m1).exec_data.wait.n_iter;
    exec_data_wait_n_reject(m1) = sim_data_record(m1).exec_data.wait.n_reject;
    exec_data_nowait_n_iter(m1) = sim_data_record(m1).exec_data.nowait.n_iter;
    exec_data_nowait_n_reject(m1) = sim_data_record(m1).exec_data.nowait.n_reject;
    
    comp_time_wait(m1) = sim_data_record(m1).compute_time.wait;
    comp_time_wait_check(m1) = sim_data_record(m1).compute_time.wait_check;
    comp_time_nowait(m1) = sim_data_record(m1).compute_time.nowait;
end

fprintf('Mean Cost diff: %f\n',mean(cost_diff))
fprintf('Max Cost diff: %f\n',max(cost_diff))

fprintf('Min Wait weight: %f\n',min(wait_weight))
fprintf('Min Move weight: %f\n',min(move_weight))
fprintf('Max Wait weight: %f\n',max(wait_weight))
fprintf('Max Move weight: %f\n',max(move_weight))

figure(1)

subplot(321)
plot(1:n_sim, cost_diff_pcent,'b')
title('Percent Cost difference: nowait - wait')
xlabel('Sim #'); ylabel('% \DeltaCost')
% legend('Cost_diff','Exposure weight')

subplot(322)
plot(min_field, cost_diff_pcent,'*')
title('Cost Diff: nowait-wait vs. Min Field Value')
xlabel('Field Min'); ylabel('% \DeltaCost')
% legend('Min Field','Wait weight, \alpha_W')

subplot(323)
semilogx(min_difft, cost_diff_pcent,'*')
title('Min Field Gradient')
xlabel('Min Gradient'); ylabel('% \DeltaCost')

subplot(324)
semilogx(max_difft, cost_diff_pcent,'*')
title('Cost difference: nowait - wait vs. Max Field Gradient')
xlabel('Max Gradient'); ylabel('% \DeltaCost')

subplot(325)
plot(avg_field, cost_diff_pcent,'*')
title('Cost difference: nowait - wait, Average Field Value')
xlabel('Avg Field'); ylabel('% \DeltaCost')

subplot(326)
plot(avg_difft, cost_diff_pcent,'*')
title('Cost difference: nowait - wait, Average Field Gradient')
xlabel('Avg Gradient'); ylabel('% \DeltaCost')

% Deeper inspection of possible parameters affecting Cost diff
figure(2)
subplot(321)
plot(1:n_sim, cost_diff_pcent,'b')
title('Cost difference: nowait - wait')
xlabel('Sim #'); ylabel('% \DeltaCost')

subplot(322)
plot(n_peaks, cost_diff_pcent,'*')
title('Number of Peaks')
xlabel('# of peaks'); ylabel('% \DeltaCost');
% legend('wait','no wait')

subplot(323)
plot(1:n_sim, path_length_wait,'b',1:n_sim, path_length_nowait,'r')
title('Wait/NoWait path length')
xlabel('Sim #');
legend('wait','no wait')

subplot(324)
plot(wait_weight, cost_diff_pcent,'*')
title('Cost diff vs. Wait weight')
xlabel('\alpha_W'); ylabel('% \DeltaCost');

subplot(326)
plot(move_weight, cost_diff_pcent,'*')
title('Cost diff vs. Move Weight')
xlabel('\alpha_M'); ylabel('% \DeltaCost');

subplot(325)
hist(cost_diff_pcent,20)
title('Percent Cost Diff')
xlabel('% \DeltaCost');

figure(3) % --------- Computation Time plots -----------------------

subplot(321)
plot(comp_time_wait, cost_diff_pcent,'*')
title('Computation Time - Waiting')
xlabel('Time (s)'); ylabel('% \DeltaCost')

subplot(323)
plot(comp_time_wait_check, cost_diff_check_pcent,'*')
title('Computation Time - Waiting with Check')
xlabel('Time (s)'); ylabel('% \DeltaCost w/ check')

subplot(325)
plot(comp_time_nowait, cost_diff_pcent,'*')
title('Computation Time - No Waiting')
xlabel('Time (s)'); ylabel('% \DeltaCost')

subplot(322)
% plot(cost_diff_pcent, comp_time_wait, '*')
hist(comp_time_wait,20)
title('Computation Time - Waiting')
xlabel('Time (s)'); %ylabel('Time (s)')

subplot(324)
% plot(cost_diff_pcent, comp_time_wait_check, '*')
hist(comp_time_wait_check,20)
title('Computation Time - Waiting with Check')
xlabel('Time (s)'); %ylabel('Time (s)')

subplot(326)
% plot(cost_diff_pcent, comp_time_nowait, '*')
hist(comp_time_nowait,20)
title('Computation Time - No Waiting')
xlabel('Time (s)');% ylabel('Time (s)')


% --- Figures to be used in paper --------------------------------------

% Weights for Wait and Move -------------------------------
figure(4)
subplot(211)
plot(wait_weight, cost_diff_pcent,'*')
set(gca,'FontSize',18)
title('Cost diff vs. Wait weight')
xlabel('\alpha_W'); ylabel('% \DeltaCost');

subplot(212)
plot(move_weight, cost_diff_pcent,'*')
set(gca,'FontSize',18)
title('Cost diff vs. Move Weight')
xlabel('\alpha_M'); ylabel('% \DeltaCost');

% Field Values --------------------------------------------
figure(5)
subplot(221)
plot(min_field, cost_diff_pcent,'*')
set(gca,'FontSize',18)
title('Min Field Value')
xlabel('Field Min'); ylabel('% \DeltaCost')

subplot(222)
plot(avg_field, cost_diff_pcent,'*')
set(gca,'FontSize',18)
title('Average Field Value')
xlabel('Avg Field'); ylabel('% \DeltaCost')
xlim([min(avg_field) max(avg_field)])

subplot(223)
h = semilogx(min_difft, cost_diff_pcent,'*');
set(gca,'FontSize',18)
title('Min Field Gradient')
xlabel('Min Gradient'); ylabel('% \DeltaCost')
xlim([min(min_difft) max(min_difft)])
% xData = get(h,'XData');
% set(gca,'Xtick',linspace(xData(1), xData(end),5))

subplot(224)
semilogx(max_difft, cost_diff_pcent,'*')
set(gca,'FontSize',18)
title('Max Field Gradient')
xlabel('Max Gradient'); ylabel('% \DeltaCost')
xlim([min(max_difft) max(max_difft)])

% Peaks -------------------------------------------------
figure(6)
plot(n_peaks, cost_diff_pcent,'*')
set(gca,'FontSize',18)
title('Number of Peaks')
xlabel('# of peaks'); ylabel('% \DeltaCost');

% Compute Time
figure(7)
subplot(221)
plot(comp_time_wait, cost_diff_pcent,'*')
set(gca,'FontSize',18)
title('Computation Time - Waiting')
xlabel('Time (s)'); ylabel('% \DeltaCost')

subplot(222)
hist(comp_time_wait,20)
set(gca,'FontSize',18)
title('Computation Time - Waiting')
xlabel('Time (s)'); %ylabel('Time (s)')

subplot(223)
plot(comp_time_nowait, cost_diff_pcent,'*')
set(gca,'FontSize',18)
title('Computation Time - No Waiting')
xlabel('Time (s)'); ylabel('% \DeltaCost')

subplot(224)
hist(comp_time_nowait,20)
set(gca,'FontSize',18)
title('Computation Time - No Waiting')
xlabel('Time (s)');% ylabel('Time (s)')
xlim([0 0.4])

% Number and Percent of Waiting cases with Full Waiting
wait_inds = find(cost_diff_pcent > 0);
fprintf('Waiting Beneficical Stats:\n')
fprintf('Percent Delta Cost > 0: %i/%i = ',size(wait_inds,2),size(cost_diff_pcent,2))
fprintf('%.2f percent\n',size(wait_inds,2)/size(cost_diff_pcent,2)*100)
wait_inds_5p = find(cost_diff_pcent > 5);
fprintf('Percent Delta Cost > 5: %i/%i = ',size(wait_inds_5p,2),size(cost_diff_pcent,2))
fprintf('%.2f percent\n',size(wait_inds_5p,2)/size(cost_diff_pcent,2)*100)
% Loose Check -----------------------------------------
figure(8)
subplot(211)
plot(comp_time_wait_check, cost_diff_check_pcent,'*')
set(gca,'FontSize',18)
title('Computation Time - Waiting with Local Check')
xlabel('Time (s)'); ylabel('% \DeltaCost w/ check')

subplot(212)
% plot(cost_diff_pcent, comp_time_wait_check, '*')
hist(comp_time_wait_check,20)
set(gca,'FontSize',18)
% title('Computation Time - Waiting with Check')
xlabel('Time (s)'); %ylabel('Time (s)')

% Number and Percent of Waiting cases with Local Test
wait_check_inds = find(cost_diff_check_pcent > 0);
fprintf('Waiting Beneficical Stats:\n')
fprintf('Percent Delta Cost > 0: %i/%i = ',size(wait_check_inds,2),size(cost_diff_check_pcent,2))
fprintf('%.2f percent\n',size(wait_check_inds,2)/size(cost_diff_check_pcent,2)*100)
wait_check_inds_5p = find(cost_diff_check_pcent > 5);
fprintf('Percent Delta Cost > 5: %i/%i = ',size(wait_check_inds_5p,2),size(cost_diff_check_pcent,2))
fprintf('%.2f percent\n',size(wait_check_inds_5p,2)/size(cost_diff_check_pcent,2)*100)

