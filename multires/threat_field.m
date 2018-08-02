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
Multiresolution path-planning with traversal costs based on a time-varying
spatial field.
	- Setup threat field
%}

function tvspf = threat_field(x, y, t, field_name, op_pointwise)

if strcmp(field_name, 'peaks')
	load tvspf_parameters_peaks.mat t_final n_peaks coeff_peaks_0 ...
		coeff_peaks_rate F_mult
elseif strcmp(field_name, 'map')
% 	load('map_q1632.mat');
    load map_q1632.mat Img;
elseif strcmp(field_name, 'poly')
	load tvspf_parameters_poly.mat t_final n_poly coeff_poly_0 ...
		coeff_poly_rate F_bias F_mult wksp
end

if op_pointwise
	tvspf	= threat_field_pt(x, y, t);
	return
end

tvspf.flag	= true;
if		strcmp(field_name, 'peaks')
	%{
		Field expression
		F(x,y,t) = 1 + F_mult*SUM_n {w_n(t) * exp( - 0.5*((x - x_n(t))/sgx_n(t))^2 - 0.5*((y - y_n(t))/sgy_n(t))^2 )}
	%}
	n_tlr	= 2;
	Fxy		= zeros(size(x,1), size(y,2), numel(t));
	dt_Fxy	= zeros(size(x,1), size(y,2), numel(t));
    dt2_Fxy	= zeros(size(x,1), size(y,2), numel(t));
	coeff_peaks_all	= zeros(5, n_peaks, numel(t));
	
	for m01 = 1:numel(t)
		Fxyt	= zeros(size(x,1), size(y,2));
		dt_Fxyt	= zeros(size(x,1), size(y,2));
        dt2_Fxyt= zeros(size(x,1), size(y,2));

		coeff_peaks_t= coeff_peaks_0 + t(m01)*coeff_peaks_rate;
		coeff_peaks_all(:, :, m01)	= coeff_peaks_t;

		for m11 = 1:n_peaks
            % Readability
            Wn = coeff_peaks_t(1, m11);
            Wn_dot = coeff_peaks_rate(1, m11);
            xn = coeff_peaks_t(2, m11);
            sigxn = coeff_peaks_t(4, m11);
            yn = coeff_peaks_t(3, m11);
            sigyn = coeff_peaks_t(5, m11);
            xn_dot = coeff_peaks_rate(2, m11);
            sigxn_dot = coeff_peaks_rate(4, m11);
            yn_dot = coeff_peaks_rate(3, m11);
            sigyn_dot = coeff_peaks_rate(5, m11);
            
			Fxytm = exp(-((x - xn).^2)./(2*sigxn^2)-((y - yn).^2)./(2*sigyn^2));
            Fxyt= Fxyt + Wn*Fxytm;
            % 1st order terms
            A11 = ((x - xn)./(sigxn^2))*xn_dot;
            A12 = ((y - yn)./(sigyn^2))*yn_dot;
            A13 = (((x - xn).^2)./(sigxn^3))*sigxn_dot;
            A14 = (((y - yn).^2)./(sigyn^3))*sigyn_dot;
            
			dt_Fxyt	= dt_Fxyt + (Wn_dot + Wn*(A11 + A12 + A13 + A14)).*Fxytm; 
            
            % Decay rates are linear, so coeff_peaks_2ndrate = 0.
            % So d^n/dt^n for n>1 coefficient rates are 0
            Wn_2dot = 0;
            xn_2dot = 0;
            yn_2dot = 0;
            sigxn_2dot = 0;
            sigyn_2dot = 0;
            % 2nd order terms
            A21 = ((x - xn)./(sigxn^2))*xn_2dot - ...
                  2*((x - xn)./(sigxn^3))*sigxn_dot + ...
                  xn_dot^2/(sigxn^2);
            A22 = ((y - yn)./(sigyn^2))*yn_2dot - ...
                  2*((y - yn)./(sigyn^3))*sigyn_dot + ...
                  yn_dot^2/(sigyn^2);
            A23 = (((x - xn).^2)./(sigxn^3))*sigxn_2dot - ...
                  3*(((x - xn).^2)./(sigxn^4))*sigxn_dot - ...
                  2*((x - xn)./(sigxn^3))*sigxn_dot;
            A24 = (((y - yn).^2)./(sigyn^3))*sigyn_2dot - ...
                  3*(((y - yn).^2)./(sigyn^4))*sigyn_dot - ...
                  2*((y - yn)./(sigyn^3))*sigyn_dot;
            
            A1 = A11 + A12 + A13 + A14; % Sum of 1st order terms
            A2 = A21 + A22 + A23 + A24; % Sum of 2nd order terms
            
            dt2_Fxyt = dt2_Fxyt + ...
                  (Wn_2dot + Wn_dot*(2*A1) + Wn*(A1.^2+A2)).*Fxytm;
		end

		Fxy(:,:,m01)	= Fxyt;
		dt_Fxy(:,:,m01) = dt_Fxyt;
        dt2_Fxy(:,:,m01) = dt2_Fxyt;
	end
	
	F_mult = 1;
	F_min	= Inf;
	F_max	= 0;
	for m01 = 1:numel(t)
		Fxy(:, :, m01)	= F_mult*Fxy(:, :, m01) + 0;					% +1 to ensure strictly positive
		
		dt_Fxy(:,:,m01) = F_mult*dt_Fxy(:, :, m01);
        dt2_Fxy(:,:,m01) = F_mult*dt2_Fxy(:, :, m01);
		
		Fxyt = Fxy(:, :, m01);
		F_min= min( F_min, min(Fxyt(:)) );
		F_max= max( F_max, max(Fxyt(:)) );
	end
		
elseif	strcmp(field_name, 'map')
	X		= double(Img.cdata);
	X_wide	= imresize(X,[(n_px + 2) (n_px + ceil(t_final*spd) + 2)], 'bilinear');
	n_colors= 256;

	X1	= X_wide(:, :, 1);		X2	= X_wide(:, :, 2);		X3	= X_wide(:, :, 3);
	cost_normalizer	= 350; %max( [max(X1(:)) max(X2(:)) max(X3(:))] );
	% col_map	= flipud(jet(n_colors));
	col_map	= (gray(n_colors));

	rgb_coeff_init =[0.17602; 0.28056; 0.87767];
	Xrgb	= rgb_coeff_init(1)*X_wide(2:(end-1), 2:(n_px + 1), 1) + ...
		rgb_coeff_init(2)*X_wide(2:(end-1), 2:(n_px + 1), 2) + ...
		rgb_coeff_init(3)*X_wide(2:(end-1), 2:(n_px + 1), 3);
	Xtrue	= Xrgb;

	
	
	rgb_coeff_fin	= [0.2990; 0.5870; 0.1140];
	rgb_coeff_t		= zeros(3, numel(t_spf));
	for m1 = 1:numel(t_spf)
		t	= t_spf(m1);
		rgb_coeff_t(:, m1)	= rgb_coeff_init + ( (t - t_init) ./ ...
			(t_final - t_init) ) * (rgb_coeff_fin - rgb_coeff_init);
	end

	Xtrue_t = zeros(2^(-jmin), 2^(-jmin), numel(t_spf));
	n_shift = 1;
	for m1 = 1:numel(t_spf)
	% 	if ~mod(m1, 3), n_shift = n_shift + 1; end
	% 	n_shift	= round( spd*t_spf(m1) );
		y0_img	= 1 + n_shift;
		yf_img	= y0_img + n_px - 1;

		Xrgb_t	= rgb_coeff_t(1, m1)*X_wide(2:(end-1), y0_img:yf_img, 1) + ...
			rgb_coeff_t(2, m1)*X_wide(2:(end-1), y0_img:yf_img, 2) + ...
			rgb_coeff_t(3, m1)*X_wide(2:(end-1), y0_img:yf_img, 3);

		Xtrue_t(:,:,m1) = Xrgb_t;
	% 	Xtrue_t(:,:,m1)	= wcodemat(Xrgb_t, n_colors);
	end

	n_tlr_order = 1;															% Degree of highest derivative considered in "Taylor" series approximation
	der_Xtrue	= zeros(n_px, n_px, n_tlr_order);

	drgb_dt		= ( 1 / (t_final - t_init) ) * (rgb_coeff_fin - rgb_coeff_init);	
	dXrgb_dt	= zeros(n_px, n_px);
	for m2 = 1:3
		dXrgb_dt = dXrgb_dt + drgb_dt(m2)*X_wide(2:(end-1), 2:(n_px+1), m2);
	end
	der_Xtrue(:, :, 1) = dXrgb_dt;
	

elseif	strcmp(field_name, 'poly')
	%{
		Field expression
		F(x,y,t) = F_bias + F_mult*SUM_m SUM_n { a_mn(t)* (x/X)^m * (y/Y)^n }
	%}
	
	F_min	= Inf;
	F_max	= 0;
	
	n_tlr	= 1;
	Fxy		= zeros(size(x,1), size(y,2), numel(t));
	dt_Fxy	= zeros(size(x,1), size(y,2), numel(t));
% 	coeff_poly_all	= zeros(n_poly, n_poly, numel(t));
	for m01 = 1:numel(t)
		Fxyt	= zeros(size(x,1), size(y,2));
		dt_Fxyt	= zeros(size(x,1), size(y,2));

		coeff_poly_t= coeff_poly_0 + t(m01)*coeff_poly_rate;
% 		coeff_poly_all(:, :, m01)	= coeff_poly_t;

		for m11 = 0:(n_poly - 1)
			for m21 = 0:(n_poly - 1)
				Fxyt	= Fxyt + coeff_poly_t(m11+1, m21+1).* ...
					(((x./wksp).^m11) .* ((y./wksp).^m21));
				
				dt_Fxyt	= dt_Fxyt + coeff_poly_rate(m11+1, m21+1).* ...
					(((x./wksp).^m11) .* ((y./wksp).^m21));
			end
		end
		
		Fxyt	= F_bias + F_mult*Fxyt;
		F_min	= min( F_min, min(Fxyt(:)) );
		F_max	= max( F_max, max(Fxyt(:)) );
		
		Fxy(:,:,m01)	= Fxyt;
		dt_Fxy(:,:,m01) = F_mult*dt_Fxyt;
	end
else
	fprintf('Unrecognized field type.\n');
	tvspf.flag	= false;
end

tvspf.field = Fxy/F_max;
if strcmp(field_name,'peaks')
    tvspf.difft	= dt_Fxy;
    tvspf.difft2 = dt2_Fxy;
else
    tvspf.difft	= dt_Fxy;
end
tvspf.min	= F_min;
tvspf.max	= F_max;
tvspf.n_tlr = n_tlr;




% tvspf.coeff_all = coeff_peaks_all;
% tvspf.coeff_rate = coeff_peaks_rate;

function tvspf_pt = threat_field_pt(x_pt, y_pt, t_pt)
	if strcmp(field_name, 'peaks')
		coeff_peaks_t_pt = coeff_peaks_0 + t_pt*coeff_peaks_rate;
		Fxyt_pt		= 0;
		dt_Fxyt_pt	= 0;
		dx_Fxyt_pt	= 0;
		dy_Fxyt_pt	= 0;
		for m011 = 1:n_peaks
			Fxytm_pt = exp( ...
				-((x_pt - coeff_peaks_t_pt(2, m011)).^2)./(2*coeff_peaks_t_pt(4, m011)^2) ...
				-((y_pt - coeff_peaks_t_pt(3, m011)).^2)./(2*coeff_peaks_t_pt(5, m011)^2) );
			Fxyt_pt= Fxyt_pt + coeff_peaks_t_pt(1, m011) * Fxytm_pt;

			dt_Fxyt_pt	= dt_Fxyt_pt + ( coeff_peaks_rate(1, m011) +  ...
				coeff_peaks_t_pt(1, m011) * (...
				+ ((x_pt - coeff_peaks_t_pt(2, m011))./(coeff_peaks_t_pt(4, m011)^2))*coeff_peaks_rate(2, m011) ...
				+ ((y_pt - coeff_peaks_t_pt(3, m011))./(coeff_peaks_t_pt(5, m011)^2))*coeff_peaks_rate(3, m011) ...
				+ (((x_pt - coeff_peaks_t_pt(2, m011)).^2)./(coeff_peaks_t_pt(4, m011)^3))*coeff_peaks_rate(4, m011) ...
				+ (((y_pt - coeff_peaks_t_pt(3, m011)).^2)./(coeff_peaks_t_pt(5, m011)^3))*coeff_peaks_rate(5, m011) ...
				) ) * Fxytm_pt;
			
			dx_Fxyt_pt	= dx_Fxyt_pt - F_mult*( (x_pt - coeff_peaks_t_pt(2, m011)) ...
				./ (coeff_peaks_t_pt(4, m011)^2) ) * coeff_peaks_t_pt(1, m011) * Fxytm_pt;
			dy_Fxyt_pt	= dy_Fxyt_pt - F_mult*( (y_pt - coeff_peaks_t_pt(3, m011)) ...
				./ (coeff_peaks_t_pt(5, m011)^2) ) * coeff_peaks_t_pt(1, m011) * Fxytm_pt;
		end

		tvspf_pt.field  = 1 + F_mult*Fxyt_pt;
		tvspf_pt.difft	= F_mult*dt_Fxyt_pt;
		tvspf_pt.diffx	= dx_Fxyt_pt;
		tvspf_pt.diffy	= dy_Fxyt_pt;
		
	elseif strcmp(field_name, 'map')
	elseif strcmp(field_name, 'poly')
		coeff_poly_t_pt = coeff_poly_0 + t_pt*coeff_poly_rate;
		Fxyt_pt		= 0;
		dt_Fxyt_pt	= 0;
		for m011 = 0:(n_poly - 1)
			for m021 = 0:(n_poly - 1)
				Fxyt_pt	= Fxyt_pt + coeff_poly_t_pt(m011+1, m021+1).* ...
					(((x_pt./wksp).^m011) .* ((y_pt./wksp).^m021));
				
				dt_Fxyt_pt	= dt_Fxyt_pt + coeff_poly_rate(m011+1, m021+1).* ...
					(((x_pt./wksp).^m011) .* ((y_pt./wksp).^m021));
			end
		end

		tvspf_pt.field  = F_bias + F_mult*Fxyt_pt;
		tvspf_pt.difft	= F_mult*dt_Fxyt_pt;
		tvspf_pt.diffx	= [];
		tvspf_pt.diffy	= [];
	end
end

end

