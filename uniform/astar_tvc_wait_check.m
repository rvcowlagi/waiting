%{
Copyright (c) 2014 Raghvendra V. Cowlagi. All rights reserved.

Copyright notice: 
=================
No part of this work may be reproduced without the written permission of
the copyright holder, except for non-profit and educational purposes under
the provisions of Title 17, USC Section 107 of the United States Copyright
Act of 1976. Reproduction of this work for commercial use is a violation of
copyright.


Disclaimer:
===========
This software program is intended for educational and research purposes.
The author and the institution with which the author is affiliated are not
liable for damages resulting the application of this program, or any
section thereof, which may be attributed to any errors that may exist in
this program.


Author information:
===================
Raghvendra V. Cowlagi, Ph.D,
Assistant Professor, Aerospace Engineering Program,
Department of Mechanical Engineering, Worcester Polytechnic Institute.
 
Higgins Laboratories, 247,
100 Institute Road, Worcester, MA 01609.
Phone: +1-508-831-6405
Email: rvcowlagi@wpi.edu
Website: http://www.wpi.edu/~rvcowlagi

The author welcomes questions, comments, suggestions for improvements, and
reports of errors in this program.

Program description:
====================
A* algorithm with time varying costs
	- Multiple goals, when any one will do
	- Terminal penalty
	- Need initial time
	- Assuming heuristic is time-invariant
	- Cycles allowed, waiting allowed: uses "expanded time-vertex pairs",
	therefore slower than the astar_tvc
	- "Discovers"  t-v adjacency on the fly
%}

function [vertex_data, idx_vt_goal, exec_data] = astar_tvc_wait_check(search_data)
%{
G				: Adjacency matrix
vertex_start	: Start node
vertex_goal	: Goal node set
heuristic		: Heuristic (N x 1 vector, where N is number of nodes)

vertex_data(n).mk	= marker, 0 = NEW, 1 = OPEN, 2 = CLOSED
vertex_data(n).d	= cost to come
vertex_data(n).b	= backpointer
vertex_data(n).t	= time step
%}

%----- Problem data in readable form
% G			= search_data.adjacency_matrix;									% A.M. of "location" vertices, not vertex-time pairs
adjacency_list = search_data.adjacency_list;
t_grid		= search_data.t_grid;
v_start		= search_data.v_start;
t_start		= search_data.t_start;
v_goal		= search_data.v_goal;											
t_goal_wnd	= search_data.t_goal_window;									
heuristic	= [search_data.heuristic; 0];									
fcn_cost	= search_data.fcn_cost;
t_step		= search_data.t_step;

do_wait_check		= strcmp(search_data.wait_check, 'strict') || strcmp(search_data.wait_check, 'loose');
wait_check_strict	= strcmp(search_data.wait_check, 'strict');

%----- Initialize struct containing search result
n_vt_pairs	= search_data.n_vertices*search_data.n_t_grid;
v_struct	= struct('mk', 0, 'd', Inf, 'b', [], 'v',[], 't', [], 'idx_t', []);
vertex_data	= repmat(v_struct, 1, n_vt_pairs + 1);							% the last one is the dummy goal

[~, idx_t_start]= ismember(t_start, t_grid);
idx_vt_start	= (idx_t_start - 1)*search_data.n_vertices + v_start;

vertex_data(idx_vt_start).mk= 1;
vertex_data(idx_vt_start).d	= 0;
vertex_data(idx_vt_start).v	= v_start;
vertex_data(idx_vt_start).t	= t_start;
vertex_data(idx_vt_start).idx_t	= idx_t_start;

n_open		= 1;
open_list	= [v_start, t_start, idx_vt_start, idx_t_start, heuristic(v_start)];
idx_vt_goal	= n_vt_pairs + 1;
goal_closed	= false;
n_iter		= 0;
n_reject	= 0;
n_vert_exp	= 1;
while (n_open ~= 0) && (~goal_closed)
	n_iter		= n_iter + 1;
	
	v_current	= open_list(1, 1);											% Get vertex from top of (sorted) open stack
	t_current	= open_list(1, 2);											% Get time of (v, t) pair from top of (sorted) open stack
	idx_vt_current	= open_list(1, 3);
	idx_t_current	= open_list(1, 4);
	vertex_data(idx_vt_current).mk = 2;										% Mark that vertex as visited
	
	n_open			= n_open - 1;
	open_list(1, :) = [];

	% Check if current time past goal window OR v_current hit dummy goal
	if (t_current >= t_goal_wnd(2)) || (v_current == 0) || (idx_t_current == search_data.n_t_grid)
		idx_nhbrs_vt = [];
	else
		if (v_current == v_goal) && (t_current >= t_goal_wnd(1)) && (t_current <= t_goal_wnd(2))
			nhbrs_v		= 0;
			idx_nhbrs_vt= idx_vt_goal;
		else
			nhbrs_v		= [v_current adjacency_list(v_current).nhbrs];		% can stop or go to neighbors			
			idx_nhbrs_vt= idx_t_current*search_data.n_vertices + nhbrs_v;
		end
    end
	
    % Evaluates the two-node local condition
	if (do_wait_check) % Waiting check is Loose or Strict
		if (wait_check_strict) && (numel(idx_nhbrs_vt)) && (numel(nhbrs_v) > 1)
			no_wait_current = true; % Skip waiting search for this node
			for m11 = adjacency_list(v_current).nhbrs
				if search_data.threat_data.cell_data(v_current + idx_t_current*search_data.n_vertices, 4) - ...
                   search_data.threat_data.cell_data(m11       + idx_t_current*search_data.n_vertices, 4) + ...
                   search_data.threat_data.cell_data(m11       + (idx_t_current + 1)*search_data.n_vertices, 4) + ...
                   search_data.wait_weight / search_data.exposure_weight < 0
                        no_wait_current = false; % if one neighbor achieves condition, do explore node for waiting
                        break;
				end
			end
		end

		if (~wait_check_strict) && (numel(idx_nhbrs_vt)) && (numel(nhbrs_v) > 1)
			no_wait_current = false;
			for m11 = adjacency_list(v_current).nhbrs
				if search_data.threat_data.cell_data(v_current + idx_t_current*search_data.n_vertices, 4) - ...
				   search_data.threat_data.cell_data(m11       + idx_t_current*search_data.n_vertices, 4) + ...
                   search_data.threat_data.cell_data(m11       + (idx_t_current + 1)*search_data.n_vertices, 4) + ...
				   search_data.wait_weight / search_data.exposure_weight >= 0
                        no_wait_current = true; % If one neighbor violates condition, don't explore node for waiting 
                        break;
				end
			end
		end
	else
		no_wait_current = false;
	end
		
	for m_nhbrs = 1:numel(idx_nhbrs_vt)										% For all (v, t) neighbors		
		v_new		= nhbrs_v(m_nhbrs); 
		t_new		= t_current + t_step;
		
		if ((v_new == v_current) && no_wait_current)
			n_reject = n_reject + 1;
			continue; % skip searching for this waiting node
		end
		
		idx_t_new	= idx_t_current + 1;
		idx_vt_new	= idx_nhbrs_vt(m_nhbrs);
		
		heuristic_current = 0;
		if idx_t_current == search_data.n_t_grid
			fprintf('At then end\n');
		end
		cost_new = fcn_cost(v_current, v_new, idx_t_current, search_data);	% Cost to go from act to new
		% Cost function must return terminal penalty as the transition cost
		% from a goal to the dummy goal.
		
		if vertex_data(idx_vt_new).mk == 0									% Unvisited
			n_vert_exp	= n_vert_exp + 1;
			vertex_data(idx_vt_new).mk	= 1;								% Mark open
			vertex_data(idx_vt_new).d	= vertex_data(idx_vt_current).d + cost_new;	% Update c2come of newly visited (v, t) pair
			vertex_data(idx_vt_new).b	= [v_current t_current idx_t_current idx_vt_current];
			vertex_data(idx_vt_new).v	= v_new;
			vertex_data(idx_vt_new).t	= t_new;
			vertex_data(idx_vt_new).idx_t	= idx_t_new;
			
			tmp_open = bin_sort(open_list(1:n_open, :), ...
				[v_new, vertex_data(idx_vt_new).t, idx_vt_new, idx_t_new, ...
				(vertex_data(idx_vt_new).d + heuristic_current)], 5);

			if numel(tmp_open) == 0
				n_open		= 0;
				open_list	= [];
			else
				n_open	= size(tmp_open, 1);
				open_list(1:n_open, :)	= tmp_open;							% Add [v_new cost] to sorted open list
			end
		elseif vertex_data(idx_vt_new).mk == 1								% Already open, update c2come if necessary
			if vertex_data(idx_vt_new).d > vertex_data(idx_vt_current).d + cost_new
				vertex_data(idx_vt_new).d	= vertex_data(idx_vt_current).d + cost_new;
				vertex_data(idx_vt_new).b	= [v_current t_current idx_t_current idx_vt_current];
				vertex_data(idx_vt_new).v	= v_new;
				vertex_data(idx_vt_new).t	= t_new;
				vertex_data(idx_vt_new).idx_t	= idx_t_new;
				
% 				[~, loc] = ismember(idx_vt_new, open_list(1:n_open, 3));
				open_list( (open_list(1:n_open, 3) == idx_vt_new), :)= [];
				n_open			= n_open - 1;
				
				tmp_open = bin_sort(open_list(1:n_open, :), ...
					[v_new, vertex_data(idx_vt_new).t, idx_vt_new, idx_t_new, ...
					(vertex_data(idx_vt_new).d + heuristic_current)], 5);
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
	
	if vertex_data(idx_vt_goal).mk	== 1, goal_closed = true; end
end
exec_data.n_iter	= n_iter;
exec_data.n_reject	= n_reject;
exec_data.n_vert_exp= n_vert_exp;
% exec_data


%**************************************************************************
function Ap = bin_sort(A, B, c)
%--------------------------------------------------------------------------
% The rows of B are inserted into A, while sorting (ascending) according to
% column c. Both A and B have the same number of columns. A is assumed to
% be sorted ascending.

[rA, cA] = size(A);
[rB, cB] = size(B);

if numel(A) == 0, Ap = B; return; end
if numel(B) == 0, return; end
% if cB ~= cA, error('A and B must have same number of columns!\n'); end

% for count = 1:rB
thisIns		= B;
thisCost	= thisIns(1, c);

if A(rA, c) <= thisCost													% Insert after last row
	Ap	= cat(1, A, thisIns);
	return;
elseif A(1, c) >= thisCost												% Insert before first row
	Ap	= cat(1, thisIns, A);
	return;
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
% 	A	= cat(1, A(1:(insRow - 1), :), thisIns, A(insRow:rA, :));
Ap	= [A(1:(insRow - 1), :); thisIns; A(insRow:rA, :)];
% end
%**************************************************************************