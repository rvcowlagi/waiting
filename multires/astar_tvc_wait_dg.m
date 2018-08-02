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

function vertex_data = astar_tvc_wait_dg(search_data)
%{
G				: Adjacency matrix
vertex_start	: Start node
vertex_goal     : Goal node set
heuristic		: Heuristic (N x 1 vector, where N is number of nodes)

vertex_data(n).mk	= marker, 0 = NEW, 1 = OPEN, 2 = CLOSED
vertex_data(n).d	= cost to come
vertex_data(n).b	= backpointer
vertex_data(n).t	= time step
%}

%----- Problem data in readable form
G			= search_data.adjacency_matrix;		% A.M. of "location" vertices, not vertex-time pairs
t_spf		= search_data.time_span;
t_step      = search_data.t_step;
v_start		= search_data.v_start;
t_start		= search_data.t_start;
vt_goal		= search_data.vt_goal;				% All vertex-time pairs (time indices from t_spf) that are goals
heuristic	= [search_data.heuristic; 0];
fcn_cost	= search_data.fcn_cost;
t_goal_wnd  = search_data.t_goal_window;
V_mr        = search_data.cell_data;

%----- Store wavelet data in search_data
search_data.n_tlr_order	= size(search_data.wvl_coeff, 1) - 1;

search_data.n_decomp	= log2(search_data.wvl_size(end, 1)/search_data.wvl_size(1, 1));	% Level of decomposition N
search_data.n_aprx_coeff= search_data.wvl_size(1, 1)^2;										% Total number of approximation coefficients
search_data.E_haar		= [1 1 1 1; 1 -1 1 -1; 1 1 -1 -1; 1 -1 -1 1];

%----- Initialize struct containing search result
n_vt_pairs	= (size(G, 1))*numel(t_spf);
v_struct	= struct('mk', 0, 'd', Inf, 'b', [], 'v',[], 't', [], 'idx_t', []);
vertex_data	= repmat(v_struct, 1, n_vt_pairs + 1);							% the last one is the dummy goal

[~, idx_t_start]= ismemberf(t_start, t_spf);
idx_vt_start	= (v_start - 1)*numel(t_spf) + idx_t_start;

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

while (n_open ~= 0) && (~goal_closed)
	n_iter		= n_iter + 1;
	
	v_current	= open_list(1, 1);					% Get vertex from top of (sorted) open stack
	t_current	= open_list(1, 2);					% Get time of (v, t) pair from top of (sorted) open stack
	idx_vt_current	= open_list(1, 3);
	idx_t_current	= open_list(1, 4);
	vertex_data(idx_vt_current).mk = 2;				% Mark that vertex as visited
	
	n_open			= n_open - 1;
	open_list(1, :) = [];

        if (idx_t_current == numel(t_spf)) || (idx_vt_current == idx_vt_goal)
            nhbrs_v		= [];
            nhbrs_t_idx = [];
            nhbrs_vt_idx= [];
		else
			if V_mr(v_current, 3) == 1	% Waiting allowed only in high-res
				nhbrs_v		= [v_current find(G(v_current, :))];			% v_current in nhbrs_v for waiting
				nhbrs_t_idx	= (idx_t_current + V_mr(v_current,3))*ones(1, numel(nhbrs_v)); % time to reach neighbor is current time +
				nhbrs_vt_idx= (nhbrs_v - 1)*numel(t_spf) + nhbrs_t_idx;                           % time to cover the distance of v_current
			else
				nhbrs_v		= find(G(v_current, :));			% v_current in nhbrs_v for waiting
				nhbrs_t_idx	= (idx_t_current + V_mr(v_current,3))*ones(1, numel(nhbrs_v)); % time to reach neighbor is current time +
				nhbrs_vt_idx= (nhbrs_v - 1)*numel(t_spf) + nhbrs_t_idx;                           % time to cover the distance of v_current
			end
			
            if any(v_current == vt_goal(:, 1)) % if v_current is one of the goals, take v_current our of neighbor list
                nhbrs_v		= nhbrs_v(2:end);  % ie. Found the goal location, but it's not time yet
                nhbrs_t_idx	= nhbrs_t_idx(2:end);
                nhbrs_vt_idx= nhbrs_vt_idx(2:end);
            end
            % If nhbrs_t_idx after time span, than no neighbors
            if any(nhbrs_t_idx > numel(t_spf))
%                 disp('reached neighbors at end of time')
                nhbrs_v		= [];
                nhbrs_t_idx = [];
                nhbrs_vt_idx= [];
            end
        end
	
	if any(vt_goal(:, 1) == v_current) && ...   % if v_current is a goal and current time inside window
        (t_current >= t_goal_wnd(1)) && (t_current <= t_goal_wnd(2))
		nhbrs_v		= [nhbrs_v (size(G, 1) + 1)]; 
		nhbrs_t_idx = [nhbrs_t_idx numel(t_spf)]; 
		nhbrs_vt_idx= [nhbrs_vt_idx idx_vt_goal];
	end

	for m_nhbrs = 1:numel(nhbrs_vt_idx)			% For all (v, t) neighbors		
		v_new		= nhbrs_v(m_nhbrs); 
		idx_t_new	= nhbrs_t_idx(m_nhbrs);
		idx_vt_new	= nhbrs_vt_idx(m_nhbrs);
		t_new		= t_spf(idx_t_new);
		
		cost_new = fcn_cost([v_current idx_t_current], [v_new idx_t_new], search_data);		% Cost to go from act to new
		% Cost function must return terminal penalty as the transition cost
		% from a goal to the dummy goal.
		
		if vertex_data(idx_vt_new).mk == 0		% Unvisited
			vertex_data(idx_vt_new).mk	= 1;	% Mark open
			vertex_data(idx_vt_new).d	= vertex_data(idx_vt_current).d + cost_new;	% Update c2come of newly visited (v, t) pair
			vertex_data(idx_vt_new).b	= [v_current t_current idx_t_current idx_vt_current];
			vertex_data(idx_vt_new).v	= v_new;
			vertex_data(idx_vt_new).t	= t_new;
			vertex_data(idx_vt_new).idx_t	= idx_t_new;
			
			tmp_open = bin_sort(open_list(1:n_open, :), ...
				[v_new, vertex_data(idx_vt_new).t, idx_vt_new, idx_t_new, ...
				(vertex_data(idx_vt_new).d + heuristic(v_new))], 5);

			if numel(tmp_open) == 0
                fprintf('A* Iteration: %i\n',n_iter)
                fprintf('temp_open == 0');
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
				
				open_list( (open_list(1:n_open, 3) == idx_vt_new), :)= [];
				n_open			= n_open - 1;
				
				tmp_open = bin_sort(open_list(1:n_open, :), ...
					[v_new, vertex_data(idx_vt_new).t, idx_vt_new, idx_t_new, ...
					(vertex_data(idx_vt_new).d + heuristic(v_new))], 5);
                if numel(tmp_open) == 0
                    fprintf('A* Iteration: %i\n',n_iter)
                    fprintf('temp_open == 0');
					n_open		= 0;
					open_list	= [];
                else
					n_open					= size(tmp_open, 1);
					open_list(1:n_open, :)	= tmp_open;		% Add [v_new cost] to sorted open list
                end
            end
		end
    end
	if vertex_data(idx_vt_goal).mk	== 2, goal_closed = true; end
end

%**************************************************************************
function A = bin_sort(A, B, c)
%--------------------------------------------------------------------------
% The rows of B are inserted into A, while sorting (ascending) according to
% column c. Both A and B have the same number of columns. A is assumed to
% be sorted ascending.

[rA, cA] = size(A);
[rB, cB] = size(B);

if numel(A) == 0, A = B; return; end
if numel(B) == 0, return; end
if cB ~= cA, error('A and B must have same number of columns!\n'); end

for count = 1:rB
	thisIns		= B(count, :);
	thisCost	= thisIns(1, c);
	
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