This folder (clean-multires) includes the main script and functions to perform multiresolution path planning and compare with the uniform resolution path. It performs the planning on the simple 3-peak field from the associated submitted paper. 

The main script to run from is "sim_discrete_mres_wait_jnl_vid.m", which also generated the submitted video of multiresolution path planning. 

To compare different resolutions, on Line 25 change jmin = -4, -5, -6, ... for increasing grid resolution as discussed in the submitted paper.