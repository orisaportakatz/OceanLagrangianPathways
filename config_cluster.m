%% solver params (PARM02)
options = odeset('RelTol',1e-3,'AbsTol',1e-3,'MaxStep',10*60*60); 
odesolver = @ode45;
%interpolant = 'spline';
interpolant = 'linear';

EOM = @EOM_polarMIT; % equation of motion dx/dt=u
K=0; % K<=0 means we don't use stochastic forcing 

ismask = 0;
isground = 0;
issteady = 0;
issharpmask = 0;
is2dim = 0;
isperiodic = 0;

ispushShore = 0;
isfreeslip = 1;

%% time params (PARM03)
% tsnapshots=; % as input to outer core
% tspan =; % as input to outer core
dt_stoch=-1;
%% grid params (PARM04)
grid_type = 'A' ; % could be 'A'/'C'
axis_dir = [1 1 1]; % each entry can be 1/-1
signOcean = 1; %inner ocean is defined with Z positive.