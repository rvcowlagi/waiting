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
	- Generate parameters for threat field 'peaks'
%}

clear all; close all; clc

t_final	= 50;
F_mult	= 1;
wksp	= 10;

coeff_peaks_0	= [...                % Fading Middle, Growing Sides 
	0.5       1.0       0.5; % weight  
   -6         0         6; % x-pos     % ** Used in FullRes exmpls 4 paper
    6         0        -6;% y-pos
	2         2         2; % x-width
	2         2         2];% y-width

coeff_peaks_f	= [...
	1.0       0.0       1.0;% weight
   -6         0         6.0;  % x-pos
    6         0        -6.0; % y-pos
	2         2         2;  % x-width
	2         2         2]; % y-width

n_peaks	= size(coeff_peaks_0,2);

coeff_peaks_rate= (coeff_peaks_f - coeff_peaks_0)/t_final;

save tvspf_parameters_peaks.mat