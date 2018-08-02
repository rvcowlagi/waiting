function cost_of_path = calc_path_cost(vertices, idx_time, search_data)

cost_of_path = 0;
for m1 = 1:numel(vertices)-1
	cost_of_path = cost_of_path + ...
		search_data.fcn_cost(vertices(m1), vertices(m1+1), idx_time(m1), search_data);
end