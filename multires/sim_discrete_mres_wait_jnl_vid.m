%{
Copyright (c) 2014. All rights reserved.

Raghvendra V. Cowlagi, Ph.D.,
Assistant Professor, Aerospace Engineering Program,
Department of Mechanical Engineering,
Worcester Polytechnic Institute.
 
Higgins Laboratories, 247,
100 Institute Road, Worcester, MA 01609.
Phone: +1-508-831-6405
Email: rvcowlagi@wpi.edu
Website: http://www.wpi.edu/~rvcowlagi

Multiresolution path-planning with traversal costs based on a time-varying 
spatial field.
	- Multiresolution in the usual way
	- Waiting allowed only in the high-resolution cells
	- Waiting costs fixed value
%}

clear variables; close all; clc;

%% Initialization
jmin	= -4;																% Coarsest possible cell has dimension 2^(-jmin)
jmax	= 0;																% Finest res cell has dimension 2^(-jmax)
n_px	= 2^(jmax-jmin);													% Total number of pixels in each row
n_decomp= -jmin;															% Coarsest level of decomposition

field_name	= 'peaks'; % Options are 'map' and 'poly'
if strcmp(field_name, 'peaks')
	load tvspf_parameters_peaks.mat
elseif strcmp(field_name, 'map')
	load('map_q1632.mat');
    t_final = 100;
    wksp = 10;
elseif strcmp(field_name, 'poly')
	load tvspf_parameters_poly.mat
end

FINE_3D_FIELD_MOVIE = false; % generate a movie of original field in 3D view
NO_ITERATIONS = false;       % Stop at first iteration of mres planner

n_pts		= 2^(-jmin);

t_init	= 0;
% Need to increase time points when gridsize/resolution increases
n_t_pts	= (t_final+1)*2^(-(jmin+4));
t_grid	= linspace(t_init, t_final, n_t_pts);

[x_wk, y_wk]= meshgrid(linspace(-wksp, wksp, n_pts), linspace(-wksp, wksp, n_pts));

tvspf		= threat_field(x_wk, y_wk, t_grid, field_name, 0);

nomcellsize	= (2*wksp) / n_px;

%% Start and Goal
p_init	= (2^(-jmin))*[0.05; 0.05];
p_goal	= (2^(-jmin))*[0.95; 0.95];

save data_tvspf_wait.mat p_init p_goal

%% Time-varying vertex data
%{
	Needed only for the computation of optimal path in the original
	(finest-resolution) graph.
%}
fprintf('Setting up finest resolution graph... \t');
tic
n_v_fine= 2^(-2*jmin);
V_fine	= zeros(n_v_fine, 4);
m		= 0;
for m1x = 1:2^(-jmin)
	for m2y = 1:2^(-jmin)
		m	= m + 1;
		V_fine(m, :)= [m1x-1 m2y-1 1 tvspf.field(m2y, m1x, 1)];
	end
end

V_fine_t	= zeros(n_v_fine, 4, numel(t_grid));
for m1 = 1:numel(t_grid)
	V_fine_t(:, :, m1) = V_fine;
	m = 0;
	for m1x = 1:2^(-jmin)
		for m2y = 1:2^(-jmin)
			m = m + 1;
			V_fine_t(m, 4, m1)	= tvspf.field(m2y, m1x, m1);
		end
	end
end

%% Finest-resolution graph
%{
	Adjacency won't change, time-varying transition costs are
	calculated in "tvcost.m"
%}
n_edges		= 0;
n_exp_edges	= n_v_fine*4;
edge_list	= zeros(n_exp_edges, 3);
for m = 1:n_v_fine
	if (m + 1 <= n_v_fine) && (mod(m, 2^(-jmin)) ~= 0)
		n_edges				= n_edges + 1;
		edge_list(n_edges, :) = [m (m + 1) 1];
		n_edges				= n_edges + 1;
		edge_list(n_edges, :) = [(m + 1) m 1];
	end

	if (m + 2^(-jmin)) <= n_v_fine
		n_edges				= n_edges + 1;
		edge_list(n_edges, :) = [m (m + 2^(-jmin)) 1];
		n_edges				= n_edges + 1;
		edge_list(n_edges, :) = [(m + 2^(-jmin)) m 1];
	end
end
G_fine = sparse(edge_list(1:n_edges,1), edge_list(1:n_edges,2), edge_list(1:n_edges,3));

p_cfine_init	= floor(p_init);
p_cfine_goal	= floor(p_goal);
[~,v_start_fine]= ismember(p_cfine_init', V_fine(:,1:2), 'rows');
[~,v_goal_fine]	= ismember(p_cfine_goal', V_fine(:,1:2), 'rows');

toc

%% Optimal path in the finest-resolution graph
search_data_fine.adjacency_matrix = G_fine;
search_data_fine.cell_data_fine_t = V_fine_t;
search_data_fine.v_start		  = v_start_fine;
search_data_fine.t_start          = t_init;	
search_data_fine.v_goal           = v_goal_fine;
search_data_fine.vt_goal          = [v_goal_fine, t_grid(end)];
search_data_fine.heuristic        = zeros(n_v_fine,1);
search_data_fine.fcn_cost         = @tvcost_mr_wait_hi_res; % cost function can handle all hi-res
search_data_fine.t_step           = ( (t_final - t_init) / (n_t_pts - 1) );
search_data_fine.mode             = 'all';
search_data_fine.const_K1         = 1e-4;
search_data_fine.time_span        = t_grid;
search_data_fine.t_goal_window    = [0, t_grid(end)];
search_data_fine.cell_data        = V_fine(:,1:3);
% -------- WAITING AND MOVEMENT COSTS -------------------------
% Set same for both Full-Res and Multi-Res
% Set wait_cost = 0.0 to enforce waiting scenarios
% Set wait_cost = Inf to enforce non-waiting scenarios
search_data_fine.wait_cost		  = 0.0;
search_data_fine.move_cost        = 0.1;
% --------------------------------------------------------------
search_data_fine.wvl_coeff        = ones(2,2);
search_data_fine.wvl_size         = ones(1,1);
% --------------------------------------------------------------
% Also need a mapping function Vmap_mr2fine that is used in 
% tvcost_mr_wait_hi_res. For fine2fine, this is just an identity function.
Vmap_fine2fine= zeros(size(V_fine,1), 1);
for m1 = 1:size(V_fine,1)
    [~, tmp2] = ismember(V_fine(m1, 1:2), V_fine(:, 1:2), 'rows');
	Vmap_fine2fine(m1) = tmp2;
end
search_data_fine.Vmap_mr2fine	 = Vmap_fine2fine;

fprintf('Searching finest resolution graph... \t');
tic_fine_start = tic;
vertex_data_fine	= astar_tvc_wait_dg(search_data_fine);
path_optimal_fine	= trace_greedy_tvc_wait_dg(search_data_fine, vertex_data_fine);
path_optimal_fine.v   = path_optimal_fine.v(1:end-1); % remove dummy goal from path
path_optimal_fine.t   = path_optimal_fine.t(1:end-1); 
tic_fine = toc(tic_fine_start);
 

%% Initialize records
tic_global_start = tic; % global TIC 

%% Multiresolution - Initialization and first iteration

visited_vertices	= zeros(n_v_fine, 1);
backpointer_fine	= zeros(n_v_fine, 1);

%----- Wavelet transform
wvlcf_all_orig			= zeros((tvspf.n_tlr + 1), 2^(-2*jmin));
[wvlcf_orig, wvlcf_sz]	= wavedec2(tvspf.field(:,:,1), n_decomp, 'db1');		% Compute wavelet decomposition of original map
wvlcf_all_orig(1, :)	= wvlcf_orig;
if strcmp(field_name,'peaks') && isfield('tvspf','difft2')
    [wvlcf_der1_orig, ~]	= wavedec2(tvspf.difft( :, :, 1), n_decomp, 'db1');
    wvlcf_all_orig(2, :)    = wvlcf_der1_orig;
    [wvlcf_der2_orig, ~]	= wavedec2(tvspf.difft2( :, :, 1), n_decomp, 'db1');
    wvlcf_all_orig(3, :)    = wvlcf_der2_orig;
else
    [wvlcf_der_orig, ~]     = wavedec2(tvspf.difft( :, :, 1), n_decomp, 'db1');
    wvlcf_all_orig(2, :)    = wvlcf_der_orig;
end
%----- MR cell decomposition
windw		= [1 1 2 2 2 2 2 2 2];
[wvlcf_all_mr, nzr_data]...
			= mrdecomposition(wvlcf_all_orig, wvlcf_sz, jmax, p_init, windw);
[G_mr, V_mr]= wvlcoeff2adjacency(nzr_data.A_nzr, wvlcf_orig, wvlcf_sz);

%----- Locate start and goal in MRCD
for j = (-n_decomp):jmax
	p_cmr_init			= floor(p_init*(2^(j)))*(2^(-j));
	[tmp, v_start_mr]	= ismember([p_cmr_init' 2^(-j)], V_mr(:,1:3), 'rows');
	if tmp, break; end
end
for j = (-n_decomp):jmax
	p_cmr_goal		= floor(p_goal*(2^(j)))*(2^(-j));
	[tmp, v_goal_mr]= ismember([p_cmr_goal' 2^(-j)], V_mr(:,1:3), 'rows');
	if tmp, break; end
end

Vmap_mr2fine= zeros(size(V_mr,1), 1);
for m1 = 1:size(V_mr,1)
% 	[~, tmp2] = ismember(V_mr(m1, 1:3), V_fine(:, 1:3), 'rows');
    [~, tmp2] = ismember(V_mr(m1, 1:2), V_fine(:, 1:2), 'rows');
	Vmap_mr2fine(m1) = tmp2;
end


%----- Setup search data struct
% search_data_mr_wait.cost_normalizer = cost_normalizer;
search_data_mr_wait.adjacency_matrix = G_mr;
search_data_mr_wait.cell_data        = V_mr;
search_data_mr_wait.v_start          = v_start_mr;
search_data_mr_wait.t_start          = t_init; %** REPEAT FOR DIFFERENT **	
search_data_mr_wait.v_goal           = v_goal_mr; %v_start_mr;
search_data_mr_wait.heuristic        = zeros(size(V_mr, 1), 1);
search_data_mr_wait.fcn_cost         = @tvcost_mr_wait_hi_res; %@tvcost_mr_wait_anywhere;
search_data_mr_wait.t_step           = ( (t_final - t_init) / (n_t_pts - 1) );
search_data_mr_wait.mode             = 'any';
search_data_mr_wait.const_K1	     = 1e-4;
search_data_mr_wait.const_K2         = 1e-4;
search_data_mr_wait.time_span        = t_grid;
search_data_mr_wait.wvl_coeff        = full(wvlcf_all_mr);
search_data_mr_wait.wvl_size         = wvlcf_sz;

search_data_mr_wait.cell_data_fine_t = V_fine_t;
search_data_mr_wait.vt_goal			 = [v_goal_mr, t_grid(end)];
search_data_mr_wait.t_goal_window    = [0, t_grid(end)];
% -------- WAITING AND MOVEMENT COSTS -------------------------
% Set above in Finest-Resolution search data
search_data_mr_wait.wait_cost		 = search_data_fine.wait_cost;
search_data_mr_wait.move_cost        = search_data_fine.move_cost;
% --------------------------------------------------------------
search_data_mr_wait.Vmap_mr2fine	 = Vmap_mr2fine;


% ---- This accounts for underestimating time/distance to goal in
% multi-resolution such as goal inside of large cell. If the goal is in a
% MR cell, then subtract the travel time of cell from time window. This
% will accomodate the extra time needed when goal is revealed in hi-res
% window. A hi-res goal uses original full time window.
search_data_mr_wait.time_span= t_grid(1:end - 2*V_mr(v_goal_mr,3)); % cell size corresponds to time idx's
search_data_mr_wait.vt_goal			= [v_goal_mr, t_grid(end) - 2*V_mr(v_goal_mr,3)*search_data_mr_wait.t_step];
search_data_mr_wait.t_goal_window   = [0, t_grid(end) - 2*V_mr(v_goal_mr,3)*search_data_mr_wait.t_step];

vertex_data_mr_wait	= astar_tvc_wait_dg(search_data_mr_wait);
toc

path_optimal_mr_wait= trace_greedy_tvc_wait_dg(search_data_mr_wait, ...
	vertex_data_mr_wait);

%% Initialize video
%-----------    Make a 3D TV Fine Field video ----------------------
if FINE_3D_FIELD_MOVIE
    video_name_tvf	= ['fine_tvspf' '_' datestr(clock,'yyyymmdd_HHMMSS.avi')];
    t_video_tvf		= 20;															% Duration of video in seconds
    video_tvf = VideoWriter(video_name_tvf);
    video_tvf.FrameRate = floor( n_t_pts/t_video_tvf );
    open(video_tvf);
    vidFig_tvf = figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.9 0.95]);
    axes_video_tvf = axes;
    minz = min(tvspf.field(:)); maxz = max(tvspf.field(:));
    for mt = 1:numel(t_grid)
        cla(axes_video_tvf);
        %axes(axes_video);
        z = tvspf.field(:, :, mt);
%         z = tvspf.difft(:, :, mt);
        surf(linspace(-wksp,wksp,n_pts),linspace(-wksp,wksp,n_pts),z)
        alpha(0.5)
        title(['Time t = ',num2str(mt)]);
        axis equal; axis tight;
        zlim([minz maxz+0.01]);
        daspect([5 5 1])
        frame_video			= getframe(vidFig_tvf);
        writeVideo(video_tvf, frame_video);
    end
    close(video_tvf);
end
% --------- Start the MR video with paths and cells ------------------
video_name	= ['mr_tvspf_wait' '_' datestr(clock,'yyyymmdd_HHMMSS.avi')];

t_video		= 20;															% Duration of video in seconds
video_global= VideoWriter(video_name);
video_global.FrameRate = floor( n_t_pts/t_video );
open(video_global);

Y = waverec2(wvlcf_all_mr(1, :), wvlcf_sz, 'db1');
z = tvspf.field(:, :, 1);
n_colors = 16;
vidFig = figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.6 0.95]);
axes_mrvideo = axes('Position',[.05 .05 .9 .9]);
colormap(gray)
colormap(flipud(colormap))

axes(axes_mrvideo)
imagesc(linspace(-wksp+0.5, wksp-0.5, n_pts), linspace(-wksp+0.5, wksp-0.5, n_pts), ...
    wcodemat(z, n_colors)); % use Y for MR or z for original field
colorbar
set(gca,'YDir','normal')
axis equal; axis tight; hold on;

axes(axes_mrvideo)
draw_cells2scale(V_mr, [], 'r', 1, 10, [], 0, nomcellsize, wksp) % Plot all Multires squares with numbers
draw_cells2scale(V_mr, [], 'b', 2, 0, path_optimal_mr_wait.v(1:end-1), 0, nomcellsize, wksp) % Draw Forward A* waiting graph

fprintf('\tIteration 0 ...\t\t\n')
fprintf('=============================================================\n');

if NO_ITERATIONS
    return
end

%% Further "incremental" planning
fprintf('Further path planning ...\t\t\n');

%----- Multi-resolution graph stuff
V_mr_current	= V_mr;
G_mr_current	= G_mr;
nzr_data_current= nzr_data;

%----- Miscellaneous
indx_p_next	= 2;
n_iter		= 0;
t_current	= 0;
is_back_step= 0;
XUp		= tvspf.field;	CfUp	= wvlcf_all_orig;
V_mr_current_draw= V_mr;
v_start_mr_next = [];
v_backptr_next = [];

path_result_fine	= [];
path_optimal_mr_current = path_optimal_mr_wait;
while V_mr_current(v_goal_mr, 3) > 1
	n_iter		= n_iter + 1;
	t_current	= t_current + search_data_mr_wait.t_step;
	
	fprintf('\tIteration %i ...\t\t\n', n_iter);
	
	% -------- Movie stuff -------
	frame_video			= getframe(vidFig); 
    writeVideo(video_global, frame_video);
	
	% -------- Step forward, note dlta = (dx, dy) --------
	if ~is_back_step
		p_current	= (V_mr_current(path_optimal_mr_current.v(indx_p_next-1), 1:2))';
		p_next		= (V_mr_current(path_optimal_mr_current.v(indx_p_next), 1:2))';
		csize_p_next= V_mr_current(path_optimal_mr_current.v(indx_p_next), 3);% Size of next cell (can force to be one)	
	else
		fprintf('backstepping... \t');
	end
	dlta		= p_next - p_current;
	
	% -------- For recording cost-to-go estimate --------
	[~,v_current_fine]	= ismember(p_current', V_fine(:,1:2), 'rows');
	[~,v_next_fine]		= ismember(p_next', V_fine(:,1:2), 'rows');
	backpointer_fine(v_next_fine)= v_current_fine;
	
	% -------- Record path and control --------
	path_result_fine	= cat(2, path_result_fine, v_current_fine);
	
	% -------- Update set of nzr detail coeffs --------
	pCells = [-n_decomp 0 0];
	for j = (-n_decomp+1):(jmax-1)											% From coarse to fine
		n = j + n_decomp + 1;
		posn	= floor((2^j)*p_current);
		pCells	= cat(1, pCells, [j posn']);
	end
	pCells = cat(1, pCells, [0 p_current']);
	nzr_data_current.pCells = pCells;	
    	
    [nzr_data_next, in_old_notin_new, in_new_notin_old]	= ...
		getNewNzrV02(nzr_data_current, wvlcf_sz, jmax, dlta, windw);

    if ~is_back_step && (numel(in_old_notin_new) == 0) && (numel(in_new_notin_old) == 0)	% In case decomposition doesn't change
		fprintf('continuing previous trajectory.\n');
        indx_p_next	= indx_p_next + 1;
		fprintf('Current location: %i\n', path_optimal_mr_current.v(indx_p_next - 1))
        fprintf('Current cell (p_current): '); disp(p_current')
        fprintf('Next cell (p_next): '); disp(p_next')
        fprintf('indx_p_next: %i\n', indx_p_next)
        disp('path_optimal_mr_current')
        disp(path_optimal_mr_current)
        % -------- Draw things --------
        cla(axes_mrvideo);
        hold on; axis equal;
        title(axes_mrvideo,['Iteration: ',num2str(n_iter),' continuing previous trajectory']) 
        z = tvspf.field(:, :, n_iter);
        % For purposes of visualizing MR graph through iterations
        p_previous	= (V_mr_current(path_optimal_mr_current.v(indx_p_next-2), 1:2))';
        [wvlcf_all_mr, nzr_data]...
			= mrdecomposition(wvlcf_all_orig, wvlcf_sz, jmax, p_previous, windw);
        Y = waverec2(wvlcf_all_mr(1, :), wvlcf_sz, 'db1');
        
        axes(axes_mrvideo)
        imagesc(linspace(-wksp+0.5, wksp-0.5, n_pts), linspace(-wksp+0.5, wksp-0.5, n_pts), ...
            wcodemat(z, n_colors)); % use Y for MR or z for original field
        colorbar
        set(gca,'YDir','normal')
        axis equal; axis tight; hold on;
        
        colormap(gray); axis equal; axis tight;
        colormap(flipud(colormap))
        
        axes(axes_mrvideo)
        draw_cells2scale(V_mr_current, [], 'r', 1, 10, [], 0, nomcellsize, wksp) % Plot all Multires squares with numbers
        draw_cells2scale(V_mr_current, [], 'b', 2, 0, path_optimal_mr_current.v(indx_p_next-1:end-1), 0, nomcellsize, wksp) % Draw Forward A* waiting graph
        
		fprintf('=============================================================\n');
        
        continue;	% Basically, don't do anything, just step fwd
        
    end % else: Decomposition changes
    indx_p_next = 2;
	
	% -------- Update graph and set of cells --------
	if is_back_step
		[wvlcf_all_mr_next, nzr_data_next]= mrdecomposition(CfUp, wvlcf_sz, jmax, p_next, windw);			% MR approximation
		[G_mr_next, V_mr_next]	= wvlcoeff2adjacency(nzr_data_next.A_nzr, wvlcf_all_mr_next, wvlcf_sz);		% V from Cf, only once; G is a matrix
		
		[CfDraw,~]	= mrDecV05(Cforig, Sz, jmax, p_next, windw);				% MR approximation
		[~, V_mr_next_draw]	= adjMat17(nzr_data_next.ANzr, CfDraw, Sz);			% V from Cf, only once; G is a matrix
	else
		[G_mr_next, V_mr_next]	= wvlcoeff2adjacency(nzr_data_current.A_nzr, wvlcf_all_orig, wvlcf_sz);
	end
	n_v_mr_next	= size(V_mr_next, 1);
	
	% -------- Locate current cell and goal in new decomposition --------
	[~, v_start_mr_next]	= ismember([p_next' csize_p_next], V_mr_next(:, 1:3), 'rows');
	for j = (-n_decomp):jmax
		p_cmr_goal		= floor(p_goal*(2^(j)))*(2^(-j));
		[tmp, v_goal_mr]= ismember([p_cmr_goal' 2^(-j)], V_mr_next(:,1:3), 'rows');
		if tmp, break; end
	end
	
	[~, v_next_mr_draw]		= ismember([p_next' csize_p_next], V_mr_next(:, 1:3), 'rows');
	for j = (-n_decomp):jmax
		p_cmr_goal_draw		= floor(p_goal*(2^(j)))*(2^(-j));
		[tmp,v_goal_mr_draw]= ismember([p_cmr_goal_draw' 2^(-j)], V_mr_next(:,1:3), 'rows');
		if tmp, break; end
	end
		
 	% -------- Redefine heuristic for changed nodes --------
	heur = zeros(n_v_mr_next,1);

	%----------------------- Run search again --------------------------
	p_cfine_backptr_next	= V_fine(backpointer_fine(v_next_fine), 1:2)';
	[~, v_backptr_next]= ismember([p_cfine_backptr_next' 1], V_mr_next(:, 1:3), 'rows');
	G_mr_next(v_start_mr_next, v_backptr_next) = 0;
	G_mr_next(v_backptr_next, v_start_mr_next) = 0;
    
    %----- Check if goal is inside high resolution window ----
    if V_mr_next(v_goal_mr, 3) == 1
       % jump outside to top of while loop, get checked again
       % and exit out to finish final path
       continue;
    end
 
    
    %----- Forwards A* from start vertex to boundary with waiting
    fprintf('Searching MRCD graph with waiting... \t\t\t');
    tic
    % Update A* with Waiting search data
    
    Vmap_mr2fine= zeros(size(V_mr_next,1), 1);
    for m1 = 1:size(V_mr_next,1)
        [~, tmp2] = ismember(V_mr_next(m1, 1:2), V_fine(:, 1:2), 'rows');
        Vmap_mr2fine(m1) = tmp2;
    end
    
    search_data_mr_wait.adjacency_matrix = G_mr_next;
    search_data_mr_wait.cell_data        = V_mr_next;
    search_data_mr_wait.heuristic        = zeros(size(V_mr_next, 1), 1);
    search_data_mr_wait.v_start          = v_start_mr_next;
    search_data_mr_wait.v_goal           = v_goal_mr;
    search_data_mr_wait.vt_goal			 = [v_goal_mr, t_grid(end)];
    search_data_mr_wait.Vmap_mr2fine	 = Vmap_mr2fine;
	search_data_mr_wait.Cf               = wvlcf_all_orig;
 	search_data_mr_wait.t_start          = t_current;	
    
    % ---- This accounts for underestimating time/distance to goal in
    % multi-resolution such as goal inside of large cell. If the goal is in a
    % MR cell, then subtract the travel time of cell from time window. This
    % will accomodate the extra time needed when goal is revealed in hi-res
    % window. A hi-res goal uses original full time window.   
    search_data_mr_wait.time_span= t_grid(1:end - 2*V_mr_next(v_goal_mr,3)); % cell size corresponds to time idx's
    search_data_mr_wait.vt_goal			= [v_goal_mr, t_grid(end) - 2*V_mr_next(v_goal_mr,3)*search_data_mr_wait.t_step];
    search_data_mr_wait.t_goal_window   = [0, t_grid(end) - 2*V_mr_next(v_goal_mr,3)*search_data_mr_wait.t_step];

    vertex_data_mr_wait	= astar_tvc_wait_dg(search_data_mr_wait);
    toc
    path_optimal_mr_wait= trace_greedy_tvc_wait_dg(search_data_mr_wait, ...
        vertex_data_mr_wait);
    disp('path_optimal_mr_wait')
    disp(path_optimal_mr_wait)

    v_data_mr_current	= vertex_data_mr_wait;
    path_optimal_mr_current= path_optimal_mr_wait;
    
    fprintf('\nCurrent location: %i\n', path_optimal_mr_current.v(indx_p_next - 1))
    fprintf('Current cell (p_current): '); disp(p_current')
    fprintf('Current cell (p_next): '); disp(p_next')
    
    % -------- Draw things --------
	cla(axes_mrvideo);
    axes(axes_mrvideo); hold on; axis equal;
    title(axes_mrvideo,['Iteration: ',num2str(n_iter)])

	z = tvspf.field(:, :, n_iter);
    % For purposes of visualizing MR graph through iterations
    [wvlcf_all_mr, nzr_data]...
			= mrdecomposition(wvlcf_all_orig, wvlcf_sz, jmax, p_current, windw);
    Y = waverec2(wvlcf_all_mr(1, :), wvlcf_sz, 'db1');
  
    axes(axes_mrvideo)
    imagesc(linspace(-wksp+0.5, wksp-0.5, n_pts), linspace(-wksp+0.5, wksp-0.5, n_pts), ...
        wcodemat(z, n_colors)); % use Y for MR or z for original field
    colorbar
    set(gca,'YDir','normal')
    axis equal; axis tight; hold on;
   
	colormap(gray); axis equal; axis tight;
    colormap(flipud(colormap))
   
    axes(axes_mrvideo)
    draw_cells2scale(V_mr_next, [], 'r', 1, 10, [], 0, nomcellsize, wksp) % Plot all Multires squares with numbers
    draw_cells2scale(V_mr_next, [], 'b', 2, 0, path_optimal_mr_current.v(indx_p_next-1:end-1), 0, nomcellsize, wksp) % Draw Forward A* waiting graph
     
	% -------- Update everything --------
	G_mr_current	= G_mr_next;
	V_mr_current	= V_mr_next;
	nzr_data_current= nzr_data_next;
fprintf('=============================================================\n');

end

%----- Forwards A* from current vertex to goal with waiting
    fprintf('Final search in Hi-res graph with waiting... \t\t\t');
    tic

    Vmap_mr2fine= zeros(size(V_mr_next,1), 1);
    for m1 = 1:size(V_mr_next,1)
        [~, tmp2] = ismember(V_mr_next(m1, 1:2), V_fine(:, 1:2), 'rows');
        Vmap_mr2fine(m1) = tmp2;
    end
    % Update Forward A* with Waiting search data
	search_data_mr_wait.Cf               = wvlcf_all_orig;
    search_data_mr_wait.adjacency_matrix = G_mr_next;
	search_data_mr_wait.cell_data        = V_mr_next;
    search_data_mr_wait.v_start          = v_start_mr_next;
    search_data_mr_wait.v_goal           = v_goal_mr;
    search_data_mr_wait.vt_goal			 = [v_goal_mr, t_grid(end)];
	search_data_mr_wait.t_start          = t_current;	
    search_data_mr_wait.Vmap_mr2fine	 = Vmap_mr2fine;
    
    % ---- Set time span and goal windows back to original time window
    % now that the goal is in hi-res window
    search_data_mr_wait.time_span= t_grid;
    search_data_mr_wait.vt_goal			= [v_goal_mr, t_grid(end)];
    search_data_mr_wait.t_goal_window   = [0, t_grid(end)];
    
    vertex_data_mr_wait	= astar_tvc_wait_dg(search_data_mr_wait);
    toc
    path_optimal_mr_wait= trace_greedy_tvc_wait_dg(search_data_mr_wait, ...
        vertex_data_mr_wait);
    % --------- End Global Timer -----------
    global_time = toc(tic_global_start);
    
    
    % --------- Draw last things -----------
    cla(axes_mrvideo);
    axes(axes_mrvideo); hold on; axis equal;
    title(axes_mrvideo,['Iteration: ',num2str(n_iter)])

	z = tvspf.field(:, :, n_iter);
 
    axes(axes_mrvideo)
    imagesc(linspace(-wksp+0.5, wksp-0.5, n_pts), linspace(-wksp+0.5, wksp-0.5, n_pts), ...
        wcodemat(z, n_colors));
    colorbar
    set(gca,'YDir','normal')
    axis equal; axis tight; hold on;
	colormap(gray); axis equal; axis tight;
    colormap(flipud(colormap))

    axes(axes_mrvideo)
    draw_cells2scale(V_mr_next, [], 'r', 1, 10, [], 0, nomcellsize, wksp) % Plot all Multires squares with numbers
    draw_cells2scale(V_mr_next, [], 'b', 2, 0, path_optimal_mr_wait.v(1:end-1), 0, nomcellsize, wksp) % Draw Forward A* waiting graph

    % -------- Movie stuff -------
	frame_video			= getframe(vidFig); 
    writeVideo(video_global, frame_video);
    
close(video_global);
diary off
%return
%%
% Add the last part of the Hi-res window trajectory to fine path

path_result_fine = [path_result_fine Vmap_mr2fine(path_optimal_mr_wait.v(1:end-1))'];

% Costs of paths
cost_fine = 0;
for m1 = 1:(numel(path_optimal_fine.v) - 1)
	cost_fine = cost_fine + tvcost_wait(path_optimal_fine.v(m1),...
                path_optimal_fine.v(m1+1), m1, search_data_fine);
end

cost_mr = 0;
for m1 = 1:(numel(path_result_fine) - 1)
	cost_mr = cost_mr + tvcost_wait(path_result_fine(m1),...
                path_result_fine(m1+1), m1, search_data_mr_wait);
end

Y = waverec2(wvlcf_all_orig(1, :), wvlcf_sz, 'db1');
figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.6 0.95]); axes_video = axes;
plot(p_init(1), p_init(2), 'wo', 'MarkerFaceColor', 'w'); hold on;
title({['Wait Cost = ',num2str(search_data_mr_wait.wait_cost),...
           ' ,Move Cost = ',num2str(search_data_mr_wait.move_cost),...
           ' ,Path time = ',num2str(numel(path_result_fine))],...
           ['Fine cost = ',num2str(cost_fine),...
           ' ,MR cost = ',num2str(cost_mr)],...
           ['Fine Calc Time = ',num2str(tic_fine),' s',...
           ' ,MR Calc Time = ',num2str(global_time),' s']},'FontSize',18);
imagesc(linspace(-wksp+0.5, wksp-0.5, n_pts), linspace(-wksp+0.5, wksp-0.5, n_pts), ...
            wcodemat(z, n_colors));
colorbar
set(gca,'YDir','normal')
axis equal; axis tight; hold on;
colormap(gray); axis equal; axis tight;
colormap(flipud(colormap))
draw_cells2scale(V_fine, [], 'r', 1, 0, path_optimal_fine.v, 'r', nomcellsize, wksp);
draw_cells2scale(V_fine, [], 'b', 2, 0, path_result_fine, 0, nomcellsize, wksp);

