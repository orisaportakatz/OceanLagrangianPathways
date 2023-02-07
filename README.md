# OceanLagrangianPathways

This batch of code uses the LPTA toolbox [1] debugged for 3D runs to calculate Lagrangian trajectories in a realistic ocean velocity field with a linear interpolator in space and time.
It is built to work with SODA3.4.2 monthly velocity field data [2], resolution 1/2°-1/2°-50 depth levels.
You can set your desired initial conditions, start month and start year, and add freeslip boundary conditions

Scheme:
1. download nc files from http://www.soda.umd.edu/ or from gw.atmos.umd.edu/.
2. run "create_velocity_as_mat_monthly.m": takes a raw SODA-monthly .nc file, extracts U,V,W velocity field, fits to mercator projection, reduces to desired region (default is North Atlantic), regrids W field to same grid as U and V, sets velocity field to 0 in the ground, and returns U,V,W as matrices with dimensions (longitude, latitude, depth, month(1-12))
3. run "fullMonthly_reg_run_cluster.m" - this runs the toolbox and eventually saves trajectories as xyzAll, where:
X = xyzAll(:,:,1);
Y = xyzAll(:,:,2);
Z = xyzAll(:,:,3);

*check that the configfile fits your requirements. set isfreeslip=1, ispushshore=0 for freeslip boundary. set interpolant='linear'.

[1]: Fredj, Erick, et al. "The particle tracking and analysis toolbox (PaTATO) for Matlab." Limnology and Oceanography: Methods 14.9 (2016): 586-599.
[2]: Carton, James A., Gennady A. Chepurin, and Ligang Chen. "SODA3: A new ocean climate reanalysis." Journal of Climate 31.17 (2018): 6967-6983.
