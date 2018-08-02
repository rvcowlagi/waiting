%{
Copyright (c) 2017. All rights reserved.

Raghvendra V. Cowlagi, Ph.D.,
Assistant Professor, Aerospace Engineering Program,
Department of Mechanical Engineering,
Worcester Polytechnic Institute.
 
Higgins Laboratories, 247,
100 Institute Road, Worcester, MA 01609.
Phone: +1-508-831-6405
Email: rvcowlagi@wpi.edu
Website: http://www.wpi.edu/~rvcowlagi

Uniform Resolution path-planning with traversal costs based on a time-varying 
spatial field.

%}

% This script compares waiting search with and without check for no-wait
% local condition
clear variables; close all; clc

dt = datestr(now,'mmmm dd, yyyy HH:MM:SS AM');
fprintf('Beginning simulation runs at: ')
disp(dt)

sim_data_record = [];

n_grid_pt	= 40;
n_vertices	= n_grid_pt^2;

speed = 1;
t_final	= 100;
wksp	= 10;

grid_sep	= 2*wksp/(n_grid_pt - 1);

t_step = grid_sep/speed; % t = dist/speed
t_max	= t_final; % from 'gen_threatfield..' or hardcode
t_grid	= 0:t_step:t_max;

N = 1;
parfor_progress(N); % Initialize progress bar
parfor n_sim_iter = 1:N
    %--------------------------------------------------------------------
    % ------ ASSIGN WEIGHTS AND FIELD PARAMETERS: MU, SIGMA, W ---------
    %--------------------------------------------------------------------
    search_data = [];
    search_data.wait_weight     = 0.1*rand;
	search_data.move_weight     = 1*rand;
	search_data.exposure_weight	= 1; %*rand;
	n_peaks		= 3 + round(37*rand);
% 	coeff_peaks_0 = [rand(1, n_peaks); ...
%                     -wksp*(-1 + 2*rand(2,n_peaks)); ...
%                     -0.25*wksp*(-1 + 2*rand(2,n_peaks))];
% 	coeff_peaks_f = [rand(1, n_peaks); ...
%                     -wksp*(-1 + 2*rand(2,n_peaks)); ...
%                     -0.25*wksp*(-1 + 2*rand(2,n_peaks))];
                
% FIXED Uniform Peak heights, VARYING location (means), FIXED widths (variances)            
%     coeff_peaks_0 = [ones(1, n_peaks); ...                  % weight
%                     -wksp*(-1 + 2*rand(2,n_peaks)); ...     % [mu_x, mu_y]'
%                     -0.25*wksp*(-1 + 2*rand(2,n_peaks))];   % [sig_x, sig_y']
%     coeff_peaks_f = [ones(1, n_peaks); ...  
%                     -wksp*(-1 + 2*rand(2,n_peaks)); ...
%                     coeff_peaks_0(4:5,:)];
%                 
% FIXED Uniform Peak heights, FIXED location (means), VARYING widths (variances)            
%     coeff_peaks_0 = [ones(1, n_peaks); ...                  % weight
%                     -wksp*(-1 + 2*rand(2,n_peaks)); ...     % [mu_x, mu_y]'
%                     -0.25*wksp*(-1 + 2*rand(2,n_peaks))];   % [sig_x, sig_y']
%     coeff_peaks_f = [ones(1, n_peaks); ...  
%                     coeff_peaks_0(2:3,:); ...
%                     -0.25*wksp*(-1 + 2*rand(2,n_peaks))];
                
% FIXED Uniform Peak heights, VARYING location (means), VARYING widths (variances)            
    coeff_peaks_0 = [ones(1, n_peaks); ...                  % weight
                    -wksp*(-1 + 2*rand(2,n_peaks)); ...     % [mu_x, mu_y]'
                    -0.25*wksp*(-1 + 2*rand(2,n_peaks))];   % [sig_x, sig_y']
    coeff_peaks_f = [ones(1, n_peaks); ...  
                    -wksp*(-1 + 2*rand(2,n_peaks)); ...
                    -0.25*wksp*(-1 + 2*rand(2,n_peaks))];
 
	coeff_peaks_rate= (coeff_peaks_f - coeff_peaks_0)/t_final;
% 	save tvspf_parameters_peaks.mat
    grid_coordinates= zeros(n_vertices, 2);
    x_grid_srch	= linspace(-wksp, wksp, n_grid_pt);
    y_grid_srch	= linspace(-wksp, wksp, n_grid_pt);
    [x_plot_grid, y_plot_grid] = meshgrid(x_grid_srch, y_grid_srch);
	% Obtain field parameters
    %---------------------------------------------------------------------
    % ------------- ASSIGN FIELD MIN AND MAX VALUES ---------------------
    %---------------------------------------------------------------------
    field_params = [];
    field_params.min_assigned = 0.1*rand;
    field_params.max_assigned = 1;
    field_params.t_final = t_final;
    field_params.n_peaks = n_peaks;
    field_params.coeff_peaks_0 = coeff_peaks_0;
    field_params.coeff_peaks_rate = coeff_peaks_rate;
	tvspf	= threat_field_args(x_plot_grid, y_plot_grid, t_grid, 'peaks', 0, field_params);

	%% Time-varying vertex data
	%{
		Needed only for the computation of optimal path in the original
		(finest-resolution) graph.
	%}
	% fprintf('Setting up finest resolution graph... \t');
	n_v_fine= n_vertices;
	V_fine	= zeros(n_v_fine, 6);
	m		= 0;
	for m2y = 1:n_grid_pt % fill horizontally bottom to top
		for m1x = 1:n_grid_pt
			m	= m + 1;
			V_fine(m, :)= [m1x-1 m2y-1 1 tvspf.field(m1x, m2y, 1) ...
				tvspf.difft(m1x, m2y, 1) tvspf.difft2(m1x, m2y, 1)];
		end
	end

	V_fine_t	= zeros(n_v_fine*numel(t_grid), 6);
	for m1 = 1:numel(t_grid)
		V_fine_t(1:n_vertices, :) = V_fine;
		m = 0;
		for m2y = 1:n_grid_pt
			for m1x = 1:n_grid_pt
				m = m + 1;
				V_fine_t( (m + (m1 - 1)*n_vertices), 4:6)	= ...
					[tvspf.field(m1x, m2y, m1) ...
					tvspf.difft(m1x, m2y, m1) tvspf.difft2(m1x, m2y, m1)];
			end
		end
	end

	% Construct fine-resolution graph
	n_edges		= 0;
	n_exp_edges	= n_vertices*4;
	edge_list	= zeros(n_exp_edges, 3);
	adjacency_list	= repmat(( struct('nhbrs', [])), n_vertices, 1);

	% Columns of grid_coordinates swapped because solutions were plotted
	% reversed over the y=x line with respect to V_fine entries
	for m = 1:n_vertices
		grid_coordinates(m, 1)	= -wksp + rem(m-1, n_grid_pt)*grid_sep;
		grid_coordinates(m, 2)	= -wksp + floor((m-1)/n_grid_pt)*grid_sep;

		if (m + 1 <= n_vertices) && (mod(m, n_grid_pt) ~= 0)
			n_edges					= n_edges + 1;
			edge_list(n_edges, :)	= [m (m + 1) 1];
			n_edges					= n_edges + 1;
			edge_list(n_edges, :)	= [(m + 1) m 1];

			adjacency_list(m).nhbrs		= [adjacency_list(m).nhbrs (m + 1)];
			adjacency_list(m + 1).nhbrs = [adjacency_list(m + 1).nhbrs m];
		end

		if (m + n_grid_pt) <= n_vertices
			n_edges					= n_edges + 1;
			edge_list(n_edges, :)	= [m (m + n_grid_pt) 1];
			n_edges					= n_edges + 1;
			edge_list(n_edges, :)	= [(m + n_grid_pt) m 1];

			adjacency_list(m).nhbrs				= [adjacency_list(m).nhbrs (m + n_grid_pt)];
			adjacency_list(m + n_grid_pt).nhbrs = [adjacency_list(m + n_grid_pt).nhbrs m];
		end

	end

	G	= sparse(edge_list(1:n_edges,1), edge_list(1:n_edges,2), edge_list(1:n_edges,3));

	%% Search paramaters bundled into one struct
    
	search_data.n_vertices				= n_vertices;
	search_data.grid_coordinates		= grid_coordinates;
	search_data.threat_data.n_peaks		= n_peaks;
	search_data.threat_data.field       = tvspf.field;
	search_data.threat_data.min         = tvspf.min;
	search_data.threat_data.max         = tvspf.max;
	search_data.threat_data.cell_data   = V_fine_t;

	search_data.t_start	= 0;
	search_data.t_grid	= t_grid;
	search_data.n_t_grid= numel(t_grid);
	search_data.t_max	= t_max;
	search_data.t_step	= t_step;

	search_data.grid_sep= grid_sep;

	search_data.adjacency_matrix= G;
	search_data.adjacency_list	= adjacency_list;
	search_data.v_start = 1; 
	search_data.v_goal	= n_grid_pt^2; 
	search_data.vt_goal = [search_data.v_goal, t_grid(end)];
	search_data.heuristic		= zeros(n_vertices, 1);
	search_data.t_goal_window	= [0 t_max];	

	search_data.fcn_cost = @calc_cost_wait_constant;

	search_data.mode = 'any';
	

	%% Search
	%----- Wait search with strict local condition check
% 	clc;
	fprintf('=========== Simulation %i ===========\n', n_sim_iter);
	
	%----- No-wait search
	tic
	[vertex_data_nowait, exec_data_nowait]	= astar_tvc(search_data);
	path_optimal_nowait	= trace_greedy_tvc(search_data.v_start, search_data.v_goal, vertex_data_nowait);
	elapsed_nowait = toc
	
	[path_added_wait, cost_path_added_wait] = ...
		calc_wait_reduction(path_optimal_nowait, search_data);
	
	search_data.wait_check	= 'loose';
	tic
	[vertex_data_wait_c, ~, exec_data_wait_c]	= astar_tvc_wait_check(search_data);
	path_optimal_wait_c	= trace_greedy_tvc_wait(search_data, vertex_data_wait_c);
	elapsed_wait_check = toc

	%----- Exhaustive wait search
	search_data.wait_check	= 'none';
	tic
	[vertex_data_wait, ~, exec_data_wait]	= astar_tvc_wait_check(search_data);
	path_optimal_wait	= trace_greedy_tvc_wait(search_data, vertex_data_wait);
	elapsed_wait = toc	

	sim_data_record(n_sim_iter).costs.wait		= vertex_data_wait(end).d;
	sim_data_record(n_sim_iter).costs.wait_check= vertex_data_wait_c(end).d;
	sim_data_record(n_sim_iter).costs.nowait	= vertex_data_nowait(search_data.v_goal).d;
	sim_data_record(n_sim_iter).costs.nowait_add= cost_path_added_wait;
	
	sim_data_record(n_sim_iter).paths.wait		= path_optimal_wait;
	sim_data_record(n_sim_iter).paths.wait_check= path_optimal_wait_c;
	sim_data_record(n_sim_iter).paths.nowait	= path_optimal_nowait;
	sim_data_record(n_sim_iter).paths.nowait_add= path_added_wait;
	
	sim_data_record(n_sim_iter).exec_data.wait		= exec_data_wait;
	sim_data_record(n_sim_iter).exec_data.wait_check= exec_data_wait_c;
	sim_data_record(n_sim_iter).exec_data.nowait	= exec_data_nowait;
    
    sim_data_record(n_sim_iter).compute_time.wait = elapsed_wait;
    sim_data_record(n_sim_iter).compute_time.wait_check = elapsed_wait_check;
    sim_data_record(n_sim_iter).compute_time.nowait = elapsed_nowait;
    
    sim_data_record(n_sim_iter).coeffs.min_assigned = field_params.min_assigned;
    sim_data_record(n_sim_iter).coeffs.max_assigned = field_params.max_assigned;
	sim_data_record(n_sim_iter).coeffs.n_peaks			= n_peaks;
	sim_data_record(n_sim_iter).coeffs.coeff_peaks_0	= coeff_peaks_0;
	sim_data_record(n_sim_iter).coeffs.coeff_peaks_f	= coeff_peaks_f;
	sim_data_record(n_sim_iter).weights.wait	= search_data.wait_weight;
	sim_data_record(n_sim_iter).weights.move	= search_data.move_weight;
	sim_data_record(n_sim_iter).weights.exposure= search_data.exposure_weight;

	% fprintf('************* No wait ************ \n')
% 	disp(path_optimal_nowait - 1)
	% disp(size(path_optimal_nowait))
	fprintf('************* No wait cost ************ \n')
	disp([vertex_data_nowait(search_data.v_goal).d ...
		calc_path_cost(path_optimal_nowait, 1:numel(path_optimal_nowait), search_data)])
% 
% 	% fprintf('************* Wait Path ************ \n')
% 	disp(path_optimal_wait.v(1:end-1) - 1)
% 	% disp(size(path_optimal_wait.v))
	fprintf('************* Wait cost ************ \n')
	disp([vertex_data_wait(end).d ...
		calc_path_cost(path_optimal_wait.v, path_optimal_wait.idx_t, search_data)])
	
	
	fprintf('************* Wait added cost ************ \n')
	disp(cost_path_added_wait)
	

	fprintf('************* Wait with loose check cost ************ \n')
	disp(vertex_data_wait_c(end).d)
	
	fprintf('\n\n')
    parfor_progress; % Count
end

%%
common_search_data = [];
common_search_data.n_vertices			   = n_vertices;

common_search_data.t_start	= 0;
common_search_data.t_grid	= t_grid;
common_search_data.n_t_grid= numel(t_grid);
common_search_data.t_max	= t_max;
common_search_data.t_step	= t_step;
common_search_data.grid_sep= grid_sep;

common_search_data.v_start = 1;
common_search_data.v_goal	= n_grid_pt^2; 
common_search_data.vt_goal = [common_search_data.v_goal, t_grid(end)];
common_search_data.heuristic		= zeros(n_vertices, 1);
common_search_data.t_goal_window	= [0 t_max];

file_name = ['tro15_fullres_wait_comparison_study02' '_' datestr(now,'yyyymmdd_HHMMSS') '.mat'];
save(file_name,'sim_data_record','common_search_data')

dt = datestr(now,'mmmm dd, yyyy HH:MM:SS AM');
fprintf('Simulation runs ended at: ')
disp(dt)

parfor_progress(0); % Clean up
