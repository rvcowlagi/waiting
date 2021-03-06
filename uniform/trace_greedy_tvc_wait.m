%{
Copyright (c) 2014.
Raghvendra V. Cowlagi, Ph.D.,
Assistant Professor, Aerospace Engineering Program,
Department of Mechanical Engineering,
Worcester Polytechnic Institute.
 
Higgins Laboratories, 247,
100 Institute Road, Worcester, MA 01609.
Phone: +1-508-831-6405
Email: rvcowlagi@wpi.edu
Website: http://www.wpi.edu/~rvcowlagi

Greedy trace of optimal path following A*: waiting involved, traces
backwards from dummy goal
%}
function path_optimal = trace_greedy_tvc_wait(search_data, vertex_data)

[~, idx_t_start]= ismember(search_data.t_start, search_data.t_grid);
idx_goal		= size(search_data.adjacency_matrix, 1) * numel(search_data.t_grid) + 1;
% idx_vt_start	= (search_data.v_start - 1)*numel(search_data.time_span) + idx_t_start;
idx_vt_start	= (idx_t_start - 1)*search_data.n_vertices + search_data.v_start;

v_current		= vertex_data(idx_goal).v;
idx_vt_current	= idx_goal;
path_optimal.v	= v_current;
path_optimal.t	= vertex_data(idx_goal).t;
path_optimal.idx_t = vertex_data(idx_goal).idx_t;

while (idx_vt_current ~= idx_vt_start)	
	v_current		= vertex_data(idx_vt_current).b(1);
	path_optimal.v	= cat(2, v_current, path_optimal.v);
	path_optimal.t	= cat(2, vertex_data(idx_vt_current).b(2), path_optimal.t);
	path_optimal.idx_t = cat(2, vertex_data(idx_vt_current).b(3), path_optimal.idx_t);
	idx_vt_current	= vertex_data(idx_vt_current).b(4);
end


