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

Transition costs for vertex-time pairs. Waiting allowed only in
full-resolution cells, therefore this function doesn't need the
wavelet-based computation of time- and cell-averaged elevation (as in
tvcost_mr).

Notation and identifier list
----------------------------
N	: level of decomposition
G	: adjacency matrix
V	: cell locations (top left vertex), dimensions, and elevations
Cf	: original coefficient array
Sz	: original coefficient size array

%}

function cost = tvcost_mr_wait_hi_res(vt1, vt2, search_data)
% vt = [vertex, time_index]
if (vt2(1) == size(search_data.adjacency_matrix,1)+1) % if v2 is the dummy goal
    cost = 0;                                         % return 0 cost transition
    return
end

%----- Retrieve wavelet coefficient data
wvlcf_all_mr= search_data.wvl_coeff;
wvlcf_sz	= search_data.wvl_size;											% size (has to do with MATLAB wavedec tool)
n_tlr_order	= search_data.n_tlr_order;

n_decomp	= search_data.n_decomp;		%log2(wvlcf_sz(end, 1)/wvlcf_sz(1, 1));		% Level of decomposition N
n_aprx_coeff= search_data.n_aprx_coeff; %wvlcf_sz(1, 1)^2;							% Total number of approximation coefficients
E_haar		= search_data.E_haar;		%[1 1 1 1; 1 -1 1 -1; 1 1 -1 -1; 1 -1 -1 1];
t_step		= search_data.t_step;

V_mr		= search_data.cell_data;

v1	= vt1(1);
v2	= vt2(1);

if (v1 == v2) && (V_mr(v1, 3) ~= 1)
	fprintf('???');
end

if V_mr(v1, 3) == 1 % in high resolution zone
	v1_fine	= search_data.Vmap_mr2fine(v1);
	v2_fine	= search_data.Vmap_mr2fine(v2);

	idx_t1	= vt1(2);
	idx_t2	= vt2(2);

	elev_c1	= search_data.cell_data_fine_t(v1_fine, 4, idx_t1);
	elev_c2	= search_data.cell_data_fine_t(v2_fine, 4, idx_t2);
else
    t01	= search_data.time_span(vt1(2));
	tf1	= t01 + (V_mr(v1, 3))*t_step;
	t02	= search_data.time_span(vt2(2));
	tf2 = t02 + (V_mr(v2, 3))*t_step;

	elev_c1	= 0;%cell_elevation(v1, 0);
	elev_c2	= 0;%cell_elevation(v2, 0);
	for t_deg = 0:n_tlr_order
		tlrcf1	= (tf1^(t_deg + 1) - t01^(t_deg + 1)) / ...
			((tf1 - t01)*factorial(t_deg + 1));		% See Eqn (16) in CDC 2014 paper
		elev_c1 = elev_c1 + tlrcf1*cell_elevation(v1, t_deg);
		
		tlrcf2	= (tf2^(t_deg + 1) - t02^(t_deg + 1)) / ...
			((tf2 - t02)*factorial(t_deg + 1));
		elev_c2 = elev_c2 + tlrcf2*cell_elevation(v2, t_deg);
	end
end

% ---------------------Cost of Waiting zone------------------------------
if (V_mr(v1, 3) == 1) && (v1 == v2)																	% Return cost of waiting at same vertex
    t01	= search_data.time_span(vt1(2));
	cost_w = search_data.wait_cost; %*elev_c1;
    tmax = numel(search_data.time_span);
    
    cost_f = V_mr(v1,3)*elev_c1;
    cost = cost_w + cost_f;
% ---------------------End of Cost of Waiting----------------------------    
else																		% Usual transition cost	
    cost_m = search_data.move_cost;
    cost_f = V_mr(v2,3)*elev_c2;
    cost = V_mr(v2,3)*cost_m + cost_f;  % scale move cost to worst cell size traversal
end

%==========================================================================
function elev_t = cell_elevation(v, elev_deg)
	%----- Cell data in (j, k, l) form
	c1_jkl	= [-log2(V_mr(v, 3)) (V_mr(v, 1:2) / V_mr(v, 3))];
	
	%----- Intensity of first cell
	k_aprx	= floor( c1_jkl(2)*( 2^(-n_decomp - c1_jkl(1)) ) );
	l_aprx	= floor( c1_jkl(3)*( 2^(-n_decomp - c1_jkl(1)) ) );
	%{
		Approximation coefficient of the "same area": for n_decomp = -jmin,
		both of these will be always 0.
	%}
% 	if elev_deg == 0
	cf_app_t= wvlcf_all_mr(elev_deg + 1, k_aprx*wvlcf_sz(1,1) + l_aprx + 1);
	cf_hor_t= wvlcf_all_mr(elev_deg + 1, n_aprx_coeff + 1 + ...
		k_aprx*wvlcf_sz(1,1) + l_aprx);
	cf_ver_t= wvlcf_all_mr(elev_deg + 1, n_aprx_coeff + 1 + ...
		k_aprx*wvlcf_sz(1,1) + l_aprx + wvlcf_sz(1,1)^2);
	cf_dgn_t= wvlcf_all_mr(elev_deg + 1, n_aprx_coeff + 1 + ...
		k_aprx*wvlcf_sz(1,1) + l_aprx + 2*wvlcf_sz(1,1)^2);

	elev_4cells	= (2^(-n_decomp))*E_haar*[cf_app_t; cf_hor_t; cf_ver_t; cf_dgn_t];

	k_prev	= k_aprx;		l_prev	= l_aprx;
	det_indx= n_aprx_coeff;
	for j_hat = (-n_decomp + 1):1:(c1_jkl(1) - 1)
		n_hat	= (j_hat + n_decomp + 1);
		det_indx= det_indx + 3*wvlcf_sz(n_hat, 1)^2;						% Starting index of detail coefficients at this level	

		k_hat	= floor( c1_jkl(2)*(2^(j_hat - c1_jkl(1))) );
		l_hat	= floor( c1_jkl(3)*(2^(j_hat - c1_jkl(1))) );

		if		(k_hat - 2*k_prev == 0) && (l_hat - 2*l_prev == 0), elev= elev_4cells(1);
		elseif	(k_hat - 2*k_prev == 0) && (l_hat - 2*l_prev == 1), elev= elev_4cells(2);
		elseif	(k_hat - 2*k_prev == 1) && (l_hat - 2*l_prev == 0), elev= elev_4cells(3);
		elseif	(k_hat - 2*k_prev == 1) && (l_hat - 2*l_prev == 1), elev= elev_4cells(4);
		else 	fprintf('????\n');
		end

		cf_hor_t= wvlcf_all_mr(elev_deg + 1, det_indx + 1 + ...
			k_hat*wvlcf_sz(n_hat + 1, 1) + l_hat);
		cf_ver_t= wvlcf_all_mr(elev_deg + 1, det_indx + 1 + ...
			k_hat*wvlcf_sz(n_hat + 1, 1) + l_hat + wvlcf_sz((n_hat + 1), 1)^2);
		cf_dgn_t= wvlcf_all_mr(elev_deg + 1, det_indx + 1 + ...
			k_hat*wvlcf_sz(n_hat + 1, 1) + l_hat + 2*wvlcf_sz((n_hat + 1), 1)^2);
		
% 		elev_4cells	= (2^j_hat)*E_haar*[elev*(2^(-j_hat)); cf_hor_t; cf_ver_t; cf_dgn_t];
		jhat_m2		= elev*(2^(-j_hat));
		elev_4cells	= (2^j_hat)*[...
			jhat_m2 + cf_hor_t + cf_ver_t + cf_dgn_t; ...
			jhat_m2 - cf_hor_t + cf_ver_t - cf_dgn_t; ...
			jhat_m2 + cf_hor_t - cf_ver_t - cf_dgn_t; ...
			jhat_m2 - cf_hor_t - cf_ver_t + cf_dgn_t];

		
		k_prev	= k_hat; 	l_prev	= l_hat;
	end
	k_hat	= c1_jkl(2);	l_hat	= c1_jkl(3);
	if		(k_hat - 2*k_prev == 0) && (l_hat - 2*l_prev == 0), elev= elev_4cells(1);
	elseif	(k_hat - 2*k_prev == 0) && (l_hat - 2*l_prev == 1), elev= elev_4cells(2);
	elseif	(k_hat - 2*k_prev == 1) && (l_hat - 2*l_prev == 0), elev= elev_4cells(3);
	elseif	(k_hat - 2*k_prev == 1) && (l_hat - 2*l_prev == 1), elev= elev_4cells(4);
	else 	fprintf('????\n');
	end
	elev_t	= elev;
	
end

end