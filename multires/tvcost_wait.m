% Created 03/13/2017, Ben Cooper
% For use in multi-resolution path with waiting conditions
% Take final fine resolution path result and output transition cost
% Check if v1 == v2 and add wait cost, else add move cost

function cost = tvcost_wait(v1, v2, t_idx, search_data)

wait_cost = search_data.wait_cost;
move_cost = search_data.move_cost;
% Cost of threat exposure at next space-time location
threat_cost = search_data.cell_data_fine_t(v2,4,t_idx+1); 

if v1 == v2
	cost = wait_cost + threat_cost;
else
	cost = move_cost + threat_cost;
end