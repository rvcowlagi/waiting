function [path_added_wait, cost_path_added_wait] = calc_wait_reduction(path_optimal_nowait, search_data)

% Checks if waiting can be added at vertex locations along a no-wait path
% to reduce its cost

n_no_wait_vertices	= numel(path_optimal_nowait);
[~, idx_t_start]	= ismember(search_data.t_start, search_data.t_grid);
path_added_wait.v	= path_optimal_nowait;
path_added_wait.t	= search_data.t_grid(idx_t_start:...
	(idx_t_start + n_no_wait_vertices - 1));
path_added_wait.idx_t  = idx_t_start: (idx_t_start + n_no_wait_vertices - 1);

m1	= 1;	% Position of vertex in no wait path
m2	= 1;	% Position of vertex with waiting
idx_t_current	= idx_t_start;
cost_current	= calc_path_cost(path_added_wait.v, path_added_wait.idx_t, search_data);
while m1 < n_no_wait_vertices - 1
% 	v_current			= path_optimal_nowait(m1);	
	wait_candidate_path.v = [path_added_wait.v(1:m2) ... 
		path_added_wait.v(m2) path_added_wait.v(m2+1:end)];
	wait_candidate_path.idx_t = [path_added_wait.idx_t(1:m2) ... 
			path_added_wait.idx_t(m2) + 1 ...
			path_added_wait.idx_t(m2 + 1:end) + 1];
		
	cost_added_1wait	= calc_path_cost(wait_candidate_path.v, ...
		wait_candidate_path.idx_t, search_data);
	
	if (  cost_added_1wait < cost_current )
		path_added_wait.v		= wait_candidate_path.v;
		path_added_wait.idx_t	= wait_candidate_path.idx_t;
		m2 = m2 + 1;		
	else
		m1 = m1 + 1;
		m2 = m2 + 1;
	end
	idx_t_current	= idx_t_current + 1;
	cost_current	= calc_path_cost(path_added_wait.v, ...
		path_added_wait.idx_t, search_data);
	
	if path_added_wait.idx_t(end) == search_data.n_t_grid - 1
		break;
	end
	
	if path_added_wait.t(end) == search_data.t_goal_window(2)
		break;
	end
end

path_added_wait.t	= search_data.t_grid(path_added_wait.idx_t);
cost_path_added_wait= calc_path_cost(path_added_wait.v, path_added_wait.idx_t, search_data);