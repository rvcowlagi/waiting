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

A* algorithm with time varying costs
	- Can handle multiple goals
	- Need initial time
	- Assuming heuristic is time-dependent
	- No cycles allowed, so "expanded time-vertex pairs" not needed
	- Cannot visit same vertex at multiple times
	- No waiting anywhere (in this search algorithm, see "astar_tvc_wait" for waiting)
%}

function [vertex_data, exec_data] = astar_tvc(search_data)
% G				: Adjacency matrix
% vertex_start	: Start node
% vertex_goal	: Goal node set
% heuristic		: Heuristic (N x 1 vector, where N is number of nodes)

% vertex_data(n).mk	= marker, 0 = NEW, 1 = OPEN, 2 = CLOSED
% vertex_data(n).d	= cost to come
% vertex_data(n).b	= backpointer
% vertex_data(n).t	= time step

% G			= search_data.adjacency_matrix;
adjacency_list = search_data.adjacency_list;
v_start		= search_data.v_start;
t_start		= search_data.t_start;
v_goal		= search_data.v_goal;
% t_goal		= search_data.t_goal;											% For multiple goals, v_goal and t_goal should be vectors of equal length
heuristic	= search_data.heuristic;
fcn_cost	= search_data.fcn_cost;
t_step		= search_data.t_step;

%----- Initialize struct containing search result
n_vertices	= size(adjacency_list, 1);
v_struct	= struct('mk', 0, 'd', Inf, 'b', [], 't', [], 'idx_t', []);
vertex_data	= repmat(v_struct, 1, n_vertices);

vertex_data(v_start).mk = 1;
vertex_data(v_start).d	= 0;
vertex_data(v_start).t	= t_start;
vertex_data(v_start).idx_t	= 1; %%***HARD CODED**


n_open		= 1;
open_list	= [v_start t_start heuristic(v_start)]; %%***HARD CODED**
goal_closed	= 0;

n_iter		= 0;
n_vert_exp	= 1;
while (n_open ~= 0) && (~goal_closed)
	n_iter		= n_iter + 1;
	v_current	= open_list(1, 1);											% Get vertex from top of (sorted) open stack
	vertex_data(v_current).mk = 2;											% Mark that vertex as visited
	
	n_open			= n_open - 1;
	open_list(1, :) = [];
	
	if vertex_data(v_current).t >= search_data.t_goal_window(2)
		nhbrs = [];
	else
		nhbrs	= adjacency_list(v_current).nhbrs;
	end
	for v_new = nhbrs														% For all neighbors		
		cost_new = fcn_cost(v_current, v_new, ...
			vertex_data(v_current).idx_t, search_data);							% Cost to go from act to new
		
		if vertex_data(v_new).mk == 0										% Unvisited
			n_vert_exp	= n_vert_exp + 1;
			vertex_data(v_new).mk	= 1;									% Mark open
			vertex_data(v_new).d	= vertex_data(v_current).d + cost_new;	% Update c2come of newly visited state
			vertex_data(v_new).b	= v_current;
			vertex_data(v_new).t	= vertex_data(v_current).t + t_step;
			vertex_data(v_new).idx_t= vertex_data(v_current).idx_t + 1;
			
			tmp_open = bin_sort(open_list(1:n_open, :), ...
				[v_new vertex_data(v_new).t ...
				(vertex_data(v_new).d + heuristic(v_new))], 3);
			%{
				OPEN list now includes (vertex, time) pairs. In this code,
				since v-cycles are not allowed, this doesn't matter.
				Vertices alone are marked open or closed, not (v, t) pairs.
				The inclusion of t in the OPEN list is for future ease of
				coding.
			%}

			if numel(tmp_open) == 0
				n_open		= 0;
				open_list	= [];
			else
				n_open	= size(tmp_open, 1);
				open_list(1:n_open, :)	= tmp_open;							% Add [v_new cost] to sorted open list
			end			
		elseif vertex_data(v_new).mk == 1									% Already open, update c2come if necessary
			if vertex_data(v_new).d > vertex_data(v_current).d + cost_new
				vertex_data(v_new).d	= vertex_data(v_current).d + cost_new;
				vertex_data(v_new).b	= v_current;
				vertex_data(v_new).t	= vertex_data(v_current).t + t_step;
				vertex_data(v_new).idx_t= vertex_data(v_current).idx_t + 1;
				
				[~, loc] = ismember(v_new, open_list(1:n_open, 1));
				open_list(loc,:)= [];
				n_open			= n_open - 1;
				
				tmp_open = bin_sort(open_list(1:n_open, :), ...
					[v_new vertex_data(v_new).t ...
					(vertex_data(v_new).d + heuristic(v_new))], 3);
				if numel(tmp_open) == 0
					n_open		= 0;
					open_list	= [];
				else
					n_open					= size(tmp_open, 1);
					open_list(1:n_open, :)	= tmp_open;						% Add [v_new cost] to sorted open list
				end
			end
		end
	end
	
	if strcmp(search_data.mode, 'any')
		goal_closed = false;
		for k = 1:numel(v_goal)
			if vertex_data(v_goal(k)).mk == 2, goal_closed = true; break; end
		end
	else
		goal_closed = true;
		for k = 1:numel(v_goal)
			if vertex_data(v_goal(k)).mk ~= 2, goal_closed = false; break; end
		end
	end
end
exec_data.n_iter	= n_iter;
exec_data.n_reject	= 0;
exec_data.n_vert_exp= n_vert_exp;
% exec_data

%**************************************************************************
function A = bin_sort(A, B, c)
%--------------------------------------------------------------------------
% The rows of B are inserted into A, while sorting (ascending) according to
% column c. Both A and B have the same number of columns. A is assumed to
% sorted ascending.

[rA, cA] = size(A);
[rB, cB] = size(B);

if numel(A) == 0, A = B; return; end
if numel(B) == 0, return; end
if cB ~= cA, error('A and B must have same number of columns!\n'); end

for count = 1:rB
	thisIns		= B(count, :);
	thisCost	= thisIns(1, c);
	
	if ismember(thisIns, A, 'rows'), 
		fprintf('This one came back!\t\t'); disp(thisIns)
		redn = redn + 1;
		continue;
	end

	if A(rA, c) <= thisCost													% Insert after last row
		A	= cat(1, A, thisIns);
		rA	= rA + 1;
		continue;
	elseif A(1, c) >= thisCost												% Insert before first row
		A	= cat(1, thisIns, A);
		rA	= rA + 1;
		continue;
	end
	
	nCand	= rA;															% Number of candidate rows in A that can have greater cost
	testRow	= 0;
	dirn	= 1;	
	while nCand > 1
		p		= floor(nCand/2);
		testRow = testRow + dirn*p;
		dirn	= 2*(A(testRow, c) < thisCost) - 1;
		nCand	= nCand - p;
	end	

	insRow = testRow + (dirn + 1)/2;										% Insert before this row in A
	A	= cat(1, A(1:(insRow - 1), :), thisIns, A(insRow:rA, :));
	rA	= rA + 1;
end
%**************************************************************************