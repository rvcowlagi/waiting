function cost = calc_cost_wait_constant(v1, v2, idx_t1, search_data)
% Spatial/temporal location dependent cost w/ constant penalty waiting
if v2 == 0 % used for dummy goal, 0 cost to transition to dummy goals
	cost = 0;
	return
end

idx_t2	= idx_t1 + 1;
if idx_t1 >= length(search_data.t_grid) % search_data.t_grid(end)
	cost = 1e9;
	return
end
cost_w = 0;
cost_m = 0;
if (v1 == v2)
    cost_w = 1; % gets multiplied by wait_weight coefficient in final cost
else
    cost_m = 1;
end


cost= search_data.wait_weight*cost_w*search_data.t_step + ...		       % waiting cost
      search_data.move_weight*cost_m*search_data.grid_sep + ...		       % movement cost
	  search_data.exposure_weight*search_data.t_step*...
	  search_data.threat_data.cell_data(v2 + ...
	  (idx_t2 - 1)*search_data.n_vertices, 4); % exposure cost


