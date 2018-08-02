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

Transition matrix from wavelet coefficients of a time-varying spatial field

Notation and identifier list
----------------------------
N	: level of decomposition
G	: adjacency matrix
V	: cell locations (top left vertex) and dimensions
Cf	: original coefficient array
Sz	: original coefficient size array

Notes
----------------------------
	Use test_adjacency_01 for testing
%}

function [G, V] = wvlcoeff2adjacency(A_nzr, Cf, Sz)

%----- Initialize
n_decomp		= log2(Sz(end, 1)/Sz(1, 1));								% Level of decomposition N
n_aprx_coeff	= Sz(1, 1)^2;												% Total number of approximation coefficients

list_add_cells	= [];														% (j, k, l) of cells to be added
list_rem_cells	= [];														% (j, k, l) of cells to be removed
%{
	In the time-varying case, the intensity must be computed during search
	as an average over a time interval specified by the search itself.
%} 

%----- Add coarse cells due to approximation coefficients
for na = 0:(n_aprx_coeff - 1)
	list_add_cells	= cat(1, list_add_cells, [(-n_decomp) ...
		floor(na/Sz(1, 1)) mod(na, Sz(1, 1))]);
end
A_nzr = sortrows(A_nzr, [1 2 3]);											% Coarser non-zero coefficients are analyzed first

%----- Apply rules 2-4 from CDC paper
for q = 1:size(A_nzr,1)
	j = A_nzr(q, 1);		k = A_nzr(q, 2);		l = A_nzr(q, 3);
	n = j + n_decomp + 1;
	
	det_indx= n_aprx_coeff;													% Index of corresponding detail coefficient
	for m = 1:(n-1)
		det_indx = det_indx + 3*Sz(m+1,1)^2;
	end	
	
	near_nzr= -n_decomp-1;
	if j ~= -n_decomp														% Find next coarser non-zero detail coefficient in the "same area"
		for jHat = (j-1):-1:(-n_decomp)
			kHat= floor( k*(2^(jHat-j)) );
			lHat= floor( l*(2^(jHat-j)) );
			
			if ismember([jHat kHat lHat], A_nzr, 'rows')
				near_nzr= jHat;
				break;
			end
		end
	end

	list_add_cells	= cat(1, list_add_cells, [(j+1)*ones(4,1) ...			% Rule 2
		 [2*k 2*k 2*k+1 2*k+1]' [2*l 2*l+1 2*l 2*l+1]']);
	
	for jHat = (near_nzr+1):(j-1)											% Rule 3
		kHat = floor( k*(2^(jHat-j)) );
		lHat = floor( l*(2^(jHat-j)) );		
		
		list_add_cells	= cat(1, list_add_cells, [(jHat+1)*ones(4,1) ...
			[2*kHat 2*kHat 2*kHat+1 2*kHat+1]' [2*lHat 2*lHat+1 2*lHat 2*lHat+1]']);	
	end	
	
	for jHat = (-n_decomp):j												% Rule 4
		kHat = floor( k*(2^(jHat-j)) );
		lHat = floor( l*(2^(jHat-j)) );
		list_rem_cells	= cat(1, list_rem_cells, [jHat kHat lHat]);
	end
end

%----- Manipulate records to get cells
list_add_cells = unique(list_add_cells, 'rows');
list_rem_cells = unique(list_rem_cells, 'rows');

[~, rem_indx] = ismember(list_rem_cells, list_add_cells, 'rows');

V_rec	= list_add_cells;
V_rec(rem_indx,:) = [];
n_cells	= size(V_rec,1);
V		= zeros(n_cells, 3);

for n = 1:n_cells
	cell_size = 2^(-V_rec(n,1));
	V(n,:) = [V_rec(n,2:3)*cell_size cell_size];
end

n_exp_edges	= 20*n_cells;	% ad hoc
list_edges	= zeros(n_exp_edges, 3);
n_edges		= 0;

for m1 = 1:n_cells
	for m2 = (m1 + 1):n_cells
		v1 = V(m1, 1:3);
		v2 = V(m2, 1:3);
		
		if v1(3) > v2(3)
			tmp = v1;
			v1	= v2;
			v2	= tmp;
		end

		indx = (v1(1:2) - v2(1:2))./v1(3);		
		if ((indx(1) <= v2(3)/v1(3)) && (indx(2) < v2(3)/v1(3))...
				&& (indx(1) >= -1) && (indx(2) > -1)) || ...
			((indx(2) <= v2(3)/v1(3)) && (indx(1) < v2(3)/v1(3))...
				&& (indx(2) >= -1) && (indx(1) > -1))						% 4-connectivity
% 			list_edges	= cat(1, list_edges, [m1 m2 1; m2 m1 1]);
			list_edges( ((n_edges+1):(n_edges+2)), :) = [m1 m2 1; m2 m1 1];
			n_edges = n_edges + 2;
		end
		
	end
end
G	= sparse(list_edges(1:n_edges, 1), list_edges(1:n_edges, 2), ...
	list_edges(1:n_edges, 3), n_cells, n_cells);


