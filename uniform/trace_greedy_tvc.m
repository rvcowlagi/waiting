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

Greedy trace of optimal path following A*: basic time-varying version no
different from regular.

%}
function path_optimal = trace_greedy_tvc(v_start, v_goal, vertex_data, v_known)

if nargin == 3
	v_current	= v_goal;
	path_optimal= v_goal;

	while (v_current ~= v_start)
		v_current	= vertex_data(v_current).b;
		path_optimal= cat(2, v_current, path_optimal);
	end
	return
end

[has_goal, indx_current] = ismember(v_goal, v_known);
if ~has_goal, error('knownNodes does not contain goal!'); end

path_optimal = v_goal;
while (indx_current ~= 1)
	indx_current	= vertex_data(indx_current).b;
	path_optimal	= cat(2, v_known(indx_current), path_optimal);
end
