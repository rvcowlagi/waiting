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

MR decomposition that forms a window around the current location. Removes
details from derivative coefficients also.
%}

function [wvlcf_all_mr, nzr_data] = mrdecomposition(wvlcf_all_orig, wvlcf_sz, jmax, p, windw)
% Cforig: original wavelet coefficients
% jmin	: coarsest cell is of dimension 2^(-jmin)
% jmax	: finest cell is of dimension 2^(-jmax)
% p		: position of agent in absolute units
% windw	: window function

% wvlcf_orig	= wvlcf_all_orig(:, :, 1);
n_tlr_order	= size(wvlcf_all_orig, 1) - 1;

%----- Initialization
n_decomp= log2(wvlcf_sz(end, 1)/wvlcf_sz(1, 1));							% Level of decomposition N
A_nzr	= [];
cell_p	= [];

%----- Add all approximation coefficients to the list
% nz_cf(1:wvlcf_sz(1, 1)^2, :) = [ones(wvlcf_sz(1,1)^2, 1) ...
% 	(1:wvlcf_sz(1, 1)^2)' (wvlcf_orig(1, 1:wvlcf_sz(1, 1)^2))'];			% List of all non-zero coefficients (1, loc, value)

nz_cf_all		= zeros( wvlcf_sz(1,1)^2, 3, n_tlr_order + 1);
for t_deg = 0:n_tlr_order
	nz_cf_all(:, :, t_deg + 1) = [ones(wvlcf_sz(1,1)^2, 1) ...
		(1:wvlcf_sz(1, 1)^2)' ...
		(wvlcf_all_orig(t_deg + 1, (1:wvlcf_sz(1, 1)^2)))'];				% List of all non-zero coefficients (1, loc, value)
end

%----- Add detail coefficients within the window
for j = (-n_decomp):(jmax-1)												% From coarse to fine
	n		= j + n_decomp + 1;
	posn	= floor((2^j)*p);
	cell_p	= cat(1, cell_p, [j posn']);
% 	fprintf('j = %i, n = %i, posn = (%i, %i)\n', j, n, posn(1), posn(2));
	for l = (posn(2) - windw(n)):(posn(2) + windw(n))
		for k = (posn(1) - windw(n)):(posn(1) + windw(n))
			ptr = wvlcf_sz(1, 1)^2;											% Index at end of appromixation coeffs
			for m = 1:(n-1)
				ptr = ptr + 3*wvlcf_sz(m+1,1)^2;
			end
			if (k>=0) && (k<wvlcf_sz(n+1,1)) && ...
					(l>=0) && (l<wvlcf_sz(n+1,1))							% if k and l within limits
				A_nzr	= cat(1, A_nzr, [j k l]);
				ptr		= ptr + k*wvlcf_sz(n+1,1) + l + 1;
% 				nz_cf	= cat(1, nz_cf, [1 ptr wvlcf_orig(ptr); ...							% horizontal
% 					1 ptr+wvlcf_sz(n+1,1)^2 wvlcf_orig(ptr+wvlcf_sz(n+1,1)^2);...			% vertical
% 					1 ptr+2*wvlcf_sz(n+1,1)^2 wvlcf_orig(ptr+2*wvlcf_sz(n+1,1)^2)]);		% diagonal
				
				nz_cf_all	= cat(1, nz_cf_all, zeros(3, 3, n_tlr_order+1));
				for t_deg = 0:n_tlr_order
					wvlcf_tlr	= wvlcf_all_orig(t_deg + 1, :);
					new_entry	= [1 ptr wvlcf_tlr(ptr); ...									% horizontal
						1 (ptr + wvlcf_sz(n+1,1)^2) (wvlcf_tlr(ptr + wvlcf_sz(n+1, 1)^2));...	% vertical
						1 (ptr + 2*wvlcf_sz(n+1,1)^2) (wvlcf_tlr(ptr + 2*wvlcf_sz(n+1, 1)^2))];	% vertical
					nz_cf_all(end-2:end, :, t_deg+1) = new_entry;
				end
				
			end
		end
	end
end

% Cf = sparse(nz_cf(:,1), nz_cf(:,2), nz_cf(:,3), 1, wvlcf_sz(end,1)^2);

for t_deg = 0:n_tlr_order
	wvlcf_tlr	= sparse(nz_cf_all(:, 1, t_deg+1), ...
		nz_cf_all(:, 2, t_deg+1), nz_cf_all(:, 3, t_deg+1), ...
		1, wvlcf_sz(end,1)^2);
	wvlcf_all_mr(t_deg+1, :) = wvlcf_tlr;
end

nzr_data.cell_p	= [cell_p; jmax floor((2^jmax)*p')];
nzr_data.A_nzr	= A_nzr;
