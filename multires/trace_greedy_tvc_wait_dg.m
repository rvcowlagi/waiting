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
function path_optimal = trace_greedy_tvc_wait_dg(search_data, vertex_data)

[~, idx_t_start]= ismemberf(search_data.t_start, search_data.time_span);

% % Regenerate idx_vt_goals from astar_tvc_wait_01
% G			= search_data.adjacency_matrix;		% A.M. of "location" vertices, not vertex-time pairs
% vt_goal		= search_data.vt_goal;				% All vertex-time pairs (time indices from t_spf) that are goals
% t_spf		= search_data.time_span;
% idx_vt_goals = [];
% for g_ind = 1:size(vt_goal,1)
%     [~,idx_gt_last] = ismember(vt_goal(g_ind,2),t_spf);
%     for tt = idx_t_start:idx_gt_last
%         idx_vt_goals = [idx_vt_goals; (vt_goal(g_ind,1)-1)*numel(t_spf)+tt];
%     end
% end
% % find which vertex data is the goal/time, search for which idx_vt_goal was
% % closed in astar_tvc_wait_01
% for k = 1:numel(idx_vt_goals)
%     if vertex_data(idx_vt_goals(k)).mk == 2
%         idx_goal = idx_vt_goals(k); 
%         break; 
%     end
% end

idx_goal		= size(search_data.cell_data, 1) * numel(search_data.time_span) + 1;
idx_vt_start	= (search_data.v_start - 1)*numel(search_data.time_span) + idx_t_start;

v_current		= vertex_data(idx_goal).v;
idx_vt_current	= idx_goal;
path_optimal.v	= v_current;
path_optimal.t	= vertex_data(idx_goal).t;

while (idx_vt_current ~= idx_vt_start)	
%    idx_vt_current % See where idx is when accessing vertex_data
   % backpointer fails
	v_current		= vertex_data(idx_vt_current).b(1);
	path_optimal.v	= cat(2, v_current, path_optimal.v);
	path_optimal.t	= cat(2, vertex_data(idx_vt_current).b(2), path_optimal.t);
	idx_vt_current	= vertex_data(idx_vt_current).b(4);
end


