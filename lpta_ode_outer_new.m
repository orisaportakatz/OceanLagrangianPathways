function [xyz,ts] = lpta_ode_outer_new(gridxyz, u, xyz0,conf_filename,TSNAPSHOT,TSPAN, Umask)
% syntax :  [xyz,ts] = lpta_ode_outer(gridxyz, u, xyz0,conf_filename,[TSNAPSHOT],[TSPAN])
% LPTA_ODE_OUTER  feeds iteratively the ode solver in lpta_ode_inner with two time snapshots between which
% the equations of motion should be integrated.
%
% Inputs
% ------
% ---- required inputs:
% - gridxyz = is a cell array of dimensions (u-components,x-components) ;
% if the grid is non-staggered, it could be (1, x-components)
% each cell in gridxyz contains a matrix with dimensions (x,y,z), denoting the positions of the velocity samples;
% 3 dimensions is just an example. one could have 1 or 2 as well.
% - u = is a cell array of dimensions (1,u-components), of the gridded velocity samples. each cell is of size (nx,ny,nz,nt)
% - xyz0 = is the tracer's initial conditions; it is a matrix (not cell array) with dimensions (n_tracers,n_dimensions)
% - conf_filename = the configuration file name (see below)
% ---- optional inputs :
% - TSNAPSHOT = times of the snapshots of velocity. If there are is only 1 number in this list  (e.g. =tsnapshot=0=), we run on a steady field.
% - TSPAN = times of dumps. in the case of a steady field run, we don't demand tspan to be within snapshots.
% both optional arguments are usually incorporated within the config file, but could be supplied through function arguments instead

%- Uind is created with velocity fields, it is of the same size as the
%velocity field for a single time-step, with 0 in ocean and 1 in ground.


% OUTPUTS:
%
% xyz = 3 dimensional matrix with dimensions [time tracer_index components], indicating the positions of the tracer in the dumping times;
% ts  = list of dump times
%
% configuration file
% ====================
% the config file is a matlab script defining several parameters. there are 3 parameter categories:
% 1. solver params
% - options     - a data structure produced by matlab's [[http://www.mathworks.com/help/matlab/ref/odeset.html][odeset]].
% - odesolver   - a function handle to your own ode solver, or one of matlab's (e.g. @ode45)
% - interpolant - a string, typically 'spline' or 'linear', to instruct our algorithm
%           how to interpolate between grid points. If the u fields include NaN's, 'spline' can't be used.
% - EOM         - a function handle, representing the equation of motion $dx/dt=EOM(x,u)$.
% In most cases, $EOM=u$. see documentation for further info
% - K           - momentum eddy diffusivity. Either a scalar, or a vector with number of elements = number of dimensions.
% 2. time params
% - tsnapshots  - times of the snapshots of velocity. If there are is only 1 number in this list  (e.g. tsnapshot=0), we run on a steady field.
% - tspan       - times of dumps. in the case of a steady field run, we don't demand tspan to be within snapshots.
% - dt_stoch    - time step for the Euler Maruyama scheme
% - tbounds     - integration bounds for the tracer equations
% 3. grid params
% - grid_type   - could be 'A'/'C' (corresponding to Arakawa's "A"- and "C"- grids)
%              or anything else (in that case, the program doesn't check the consistency
%              of your initial conditions with the grid definition)
% - axis_dir    - a vector of 1's and -1's, indicating the arrangement of each axis in the matrix
%                  (from small to big indices), with respect to "normal" right handed system,
%                  where north,east and up are positive. for example, MITgcm's axis_dir is [1 1 -1]
%                  because the first vertical layer is the highest and the last one is the lowest,
%                  opposite to the right handed system.
%----------------------------------------------------------------------------------------------------
% Author: Avi Gozolchiani
% $Log: lpta_ode_outer.m,v $
% Revision 1.2  2015/03/09 15:58:33  avigoz
% adding time range (tbounds) and copyright headers
%
% Revision 1.7  2015/03/09 15:47:18  avigoz
% Summary: including tbounds - beginning and end time of the particles
%
%----------------------------------------------------------------------------------------------------
run(conf_filename);
%% defaults
if ~exist( 'odesolver', 'var' )
    odesolver = @ode45;
end

%Ori - default sign of particles inside ocean is +1:
if ~exist('signOcean', 'var')
    signOcean = 1;
end

if ~exist( 'options', 'var' )
    options = [];
end
if ~exist( 'tsnapshots', 'var' )
    tsnapshots = TSNAPSHOT;
end
if ~exist( 'tspan', 'var' )
    tspan = TSPAN;
end
nsnaps=length(tsnapshots);
n_tracers=size(xyz0,1);
flag_steady=false;
if(nsnaps==1)
    flag_steady=true;
end % if(nsnaps==1)
%% some input checking :
n_components=length(u);
if(n_components>3)
    error('more than 3 dimensions is not currently covered');
end
if(~iscell(gridxyz))
    error('gridxyz should be a cell array');
end
% if it's a non staggerred grid, duplicate the mesh:
% (griddedInterpolant doesn't multiply the mesh internally,
% so there's no memory penalty)
if(any(size(gridxyz)==1))
   gridxyz=repmat(gridxyz(:)',[length(gridxyz),1]);
end % if(any(size(gridxyz)==1))

if(~all(size(gridxyz)==n_components | size(gridxyz)==1))
    error('inconsistent number of components');
end % if(~all(size(gridxyz)==n_components))
if(length(axis_dir)~=n_components)
    error('inconsistent number of axes and components');
end % if(length(axis_dir)~=n_components)
sz_grd=cell(size(gridxyz));
sz_u=cell(size(u));
for i_component_u=1:n_components
    for i_component_r=1:n_components
        sz_grd{i_component_u,i_component_r}=size(gridxyz{i_component_u,i_component_r});
        sz_u{i_component_u}=size(u{i_component_u});
        if(length(sz_u{i_component_u})==max(2,n_components)) % if the time dimension is singleton
            % in 1 dimension, sz_u has a singleton dimension anyway, so we never
            % satisfy length(sz_u{i_component_u})==1
            sz_u{i_component_u}=[sz_u{i_component_u} 1];
        end
        if(sz_u{i_component_u}(end)~=nsnaps)
            error('last dimension of velocities matrix should be equal to the number of snapshots');
        end % if(size(u,end)~=nsnaps)
        if(~all(sz_grd{i_component_u,i_component_r}==sz_u{i_component_u}(1:end-1)))
            disp([sz_grd{i_component_u,i_component_r}]); disp([sz_u{i_component_u}]);
            error('domain and velocity dimensions do not match')
        end % if(~all(sz1==sz2) | ...
    end % for i_component_r=1:n_components
end % for i_component_u=1:n_components

if(size(xyz0,2)~=n_components)
    error('the size of the initial conditions in xyz0 is wrong');
end % if(size(xyz0,2)~=n_components)
ranges_tracers=[min(xyz0);max(xyz0)];

% get the resolution at the ends and the last grid points at those ends
sep={};
ends={};
dims=1:n_components;
for i_dim=1:n_components
    sep{i_dim}=resolution_ends(gridxyz(:,dims(i_dim)),dims(i_dim),dims([1:i_dim-1 i_dim+1:end]));
    ends{i_dim}=ranges(gridxyz(:,dims(i_dim)),dims(i_dim),dims([1:i_dim-1 i_dim+1:end]));
    if(diff(ends{i_dim})<0)
        ends{i_dim}=fliplr(ends{i_dim});
        sep{i_dim}=fliplr(sep{i_dim});
    end
end

% on a non-staggered grid, the end points are actually the middle of a cell face, and we need to add one more half cell.
% on a staggered grid, we don't need to do it
range_check_flag=false;
add_half_cell=0;
switch upper(grid_type)
    case 'A'
        range_check_flag=true;
        add_half_cell=cellfun(@(x)[-0.5 0.5].*x,sep,'UniformOutput',false);
    case 'C'
        range_check_flag=true;
        add_half_cell=cellfun(@(x)[-0 0],sep,'UniformOutput',false);
end % switch upper(grid_type)

% now that we collected the extreme limits of the grid in all dimensions, we can check if the init conditions are actually within the grid. init conditions may still be on land points, though.
% that's the user's responsibility to initialize in wet point, otherwise the velocity field is just zero, and the tracer does not move,
% which is perfectly legitimate from a tracer code point of view.
if range_check_flag % if you use a kind of grid we don't know, we can't check the validity of your ranges. tracing might still be fine, though...
    ranges_grid=cellfun(@(i)(ends{i}+add_half_cell{i}),num2cell(dims),'UniformOutput',false)';
    ranges_grid=cell2mat(ranges_grid)';
    if(any(ranges_tracers(1,:)<ranges_grid(1,:)) | ...
            any(ranges_tracers(2,:)>ranges_grid(2,:)))
        disp('ranges_tracers:')
        disp(ranges_tracers);
        disp('ranges_grid:')
        disp(ranges_grid);
        error('tracers are not within the grid bounds');
    end % if(any(ranges_tracers(1,:)<ranges_grid(1,:)) | ...
end % if range_check_flag

if(~issorted(tspan))
    error('tspan is not sorted');
end % if(~issorted(tspan))

% end of input checking

%% standardize
% for non-conventional grid ordering (north/east/up are positive),
% flip the necessary dimensions in order to make it a conventional right handed system
for i_axis=1:n_components
    if(axis_dir(i_axis)<0)
        gridxyz=cellfun(@(x)(flipdim(x,i_axis)),gridxyz,'UniformOutput',false);
        u      =cellfun(@(x)(flipdim(x,i_axis)),u      ,'UniformOutput',false);
    end % if(axis_dir(i_axis)<0)
end % for i_axis=1:n_components
% K is defined for each axis
if(K(1)>0)
    if(length(K)==1)
        K=repmat(K,n_components,1);
    end
end

%% back - trajectory option
if exist( 'back_flag', 'var' )
    if back_flag
        u =cellfun(@(x)(flipdim(x,n_components+1)),u ,'UniformOutput',false);
        u =cellfun(@(x) -1.*x, u, 'UniformOutput', false);
        tspan = flipdim(tsnapshots(end)+tsnapshots(1)-tspan,2);
        if exist( 'tbounds', 'var' )
            tbounds = flipdim(tsnapshots(end)+tsnapshots(1)-tbounds,2);
        end %exist( 'tbounds', 'var' )
    end %back_flag
end %exist( 'back_flag', 'var' )


% build indices for slicing the elements in a velocity cell array
idx.type='()';                  % indices structure
idx.subs=repmat({':'},1,n_components+1);
% sliceu=@(x,idx_iters)subsref(x,idx_iters);
if exist( 'tbounds', 'var' )
    from_snap=(tsnapshots>tbounds(1));
    from_snap(1:end-1)=from_snap(1:end-1) | from_snap(2:end);
    to_snap=(tsnapshots<tbounds(2));
    to_snap(2:end)=to_snap(2:end) | to_snap(1:end-1);
    relevant_snaps=(from_snap & to_snap);
    idx.subs{n_components+1}=relevant_snaps; % relevant time slice
    u=cellfun(@(x)subsref(x,idx),u,'UniformOutput',false);
    tsnapshots=tsnapshots(relevant_snaps);
    nsnaps=sum(relevant_snaps);
    % linearly interpolate in time the first and last snaps to the beginning and end point of integration
    % u(t)=(u1-u0)/dt*(t-t0)+u0 =  u1/dt*(t-t0)+u0*(1-(t-t0)/dt) = a1*u1+a0*u0
    % where a1=(t-t0)/dt , a0=(1-(t-t0)/dt)
    interp_indices=[1 nsnaps-1];
    for i_interp=1:length(interp_indices)
        a1=(tbounds(i_interp)-tsnapshots(interp_indices(i_interp)))/diff(tsnapshots(interp_indices(i_interp)+[0 1]));
        a0=1-a1;
        idx.subs{n_components+1}=interp_indices(i_interp); % relevant time slice
        tmp_u=cellfun(@(x)subsref(x,idx)*a0,u,'UniformOutput',false);
        idx.subs{n_components+1}=interp_indices(i_interp)+1; % relevant time slice
        tmp_u=cellfun(@(x,y)(y+subsref(x,idx)*a1),u,tmp_u,'UniformOutput',false);
        idx.subs{n_components+1}=interp_indices(i_interp)+i_interp-1; % relevant time slice
        u=cellfun(@(x,y)subsasgn(x,idx,y),u,tmp_u,'UniformOutput',false);
    end
    tsnapshots(1)=tbounds(1);
    tsnapshots(nsnaps)=tbounds(2);
end % if ~exist( 'tbounds', 'var' )
Org_tsnapshots=tsnapshots;
%% call the inner iterations :

cur_xyz0=xyz0;
xyz=reshape(xyz0,[1 size(xyz0)]);
ts=tsnapshots(1);
if(flag_steady)
    disp('here')
    [xyz,ts] = lpta_ode_inner_new(gridxyz, u, [0 0],tspan, xyz0, interpolant,EOM ,K,dt_stoch,options, odesolver, signOcean ,Umask);
end

for i_iter=1:nsnaps-1 % every iteration interpolates between two snapshots, so n_iterations=nsnaps-1
    t2pts=tsnapshots([0 1]+i_iter);
    cur_tspan=tspan(tspan>=t2pts(1) & tspan<=t2pts(2));
    cur_tspan2=unique([t2pts(1) cur_tspan(:)' t2pts(2)]); % add the snapshot points in case they don't exist in the users's dump requests
    idx.subs{n_components+1}=[i_iter i_iter+1]; % relevant time slice
    [tmpxyz,tmpts] = lpta_ode_inner_new(gridxyz, cellfun(@(x)subsref(x,idx),u,'UniformOutput',false), t2pts,cur_tspan2, cur_xyz0, interpolant,EOM ,K,dt_stoch,options, odesolver, signOcean ,Umask);
    %%%%%%
    % 0/360 boundery (Israela 5.2019)- checking
    if range_check_flag
        n_tracer=length(tmpxyz(1,:,1));
        for i=1:n_tracer
            index_dir_E=find(tmpxyz(:,i,1)>=ranges_grid(2,1));
            index_dir_W=find(tmpxyz(:,i,1)<=ranges_grid(1,1));
            tmpxyz(index_dir_E,i,1)=tmpxyz(index_dir_E,i,1)-360;
            tmpxyz(index_dir_W,i,1)=tmpxyz(index_dir_W,i,1)+360;
        end
    end
    %%%%%%%
    cur_xyz0=squeeze(tmpxyz(end,:,:));
    % if the user didn't want the dumps in the snapshot points, remove them
    [~,subset_t,~]=intersect(tmpts,cur_tspan);
    tmpxyz=tmpxyz(subset_t,:,:);
    tmpts=tmpts(subset_t);
    % check ranges
    
    ranges_tracers=[min(tmpxyz,[],1);max(tmpxyz,[],1)];
    if range_check_flag % if you use a kind of grid we don't know, we can't check the validity of your ranges. tracing might still be fine, though...
        Ioutrange=false(1,n_tracers);%indicator function
        for i_component=1:n_components
            Ioutrange=Ioutrange | ...
                (ranges_tracers(1,:,i_component)<ranges_grid(1,i_component));
            Ioutrange=Ioutrange | ...
                (ranges_tracers(2,:,i_component)>ranges_grid(2,i_component));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %HERE ALL PARTICLES THAT LEFT THE RANGE ARE DEALT WITH.
        %ORIGINAL SOLUTION: SET VALUE TO NAN:
        %     tmpxyz(:,Ioutrange,:)=nan;
        %NEW SOLUTION (27.6.19): SET VALUE TO SHORE
        %tmpxyz(:,Ioutrange,:) =
        tmpxyz(:,Ioutrange,1)= -100;
        tmpxyz(:,Ioutrange,2)= 50;
        tmpxyz(:,Ioutrange,3)= 100;
        
    end % if range_check_flag
    cur_xyz0=squeeze(tmpxyz(end,:,:));
    
    % concatenate the results of different snapshots
    xyz=cat(1,xyz,tmpxyz);
    ts=cat(1,ts,tmpts);
end % for i_iter=1:nsnaps-1
[ts,its,~]=unique(ts);
xyz=xyz(its,:,:);

%% arrange ts to be backward if back-trajectories
if exist( 'back_flag', 'var' )
    if back_flag
        ts=Org_tsnapshots(end)+Org_tsnapshots(1)-ts;
    end %back_flag
end %exist( 'back_flag', 'var' )

%%%%%%----------------------Subfunctions--------------------%%%%%%%
function sep=resolution_ends(gridxyz,cur_dir_index,cross_section_dir_index)
% sufficient but not necessary condition to avoiding leakage:
% for each y (without loss of generality) take the maximal x-resolution dx in both (east and west)
% ends of the grid. In cases of a staggered horizontal grid, maximize these
% two values among the different cell definitions (according to the
% different staggered variables). Maximize among all y- cross sections as well
% gridxyz - is a column matrix of cells corresponding to the axis we are interested in .
% for example, for y this might be gridxyz(:,2) of particle_track_ode_outer
% cur_dir_index - direction of the resolution we want to calculate
% cross_section_dir_index - orthogonal to cur_dir_index.
%
sep=cellfun(@(x)(maxndim(abs(diff(x,1,cur_dir_index)),cross_section_dir_index)),gridxyz,'UniformOutput',false); % max among all y- cross sections (for example)
sep=cellfun(@(x)([x(1) x(end)]),sep,'UniformOutput',false); % two ends
sep=cell2mat(sep);
sep=max(sep,[],1); % max among staggered variables

function ends=ranges(gridxyz,cur_dir_index,cross_section_dir_index)
% sufficient but not necessary condition to avoiding leakage of tracers:
% for each y (without loss of generality) take the maximal x-resolution dx in both (east and west)
% ends of the grid. In cases of a staggered horizontal grid, maximize these
% two values among the different cell definitions (according to the
% different staggered variables). Maximize among all y- cross sections as well
% gridxyz - is a column matrix of cells corresponding to the axis we are interested in .
% for example, for y this might be gridxyz(:,2) of particle_track_ode_outer
% cur_dir_index - direction of the resolution we want to calculate
% cross_section_dir_index - orthogonal to cur_dir_index.

% maxndim is maximum over multiple dimensions, -maxndim(-x) is
% the minimum ...
end1=cellfun(@(x)(-maxndim(-x,cross_section_dir_index)),gridxyz,'UniformOutput',false); % two ends
end1=cellfun(@(x)(x(1)),end1,'UniformOutput',false); % two ends
end2=cellfun(@(x)(maxndim(x,cross_section_dir_index)),gridxyz,'UniformOutput',false); % two ends
end2=cellfun(@(x)(x(end)),end2,'UniformOutput',false); % two ends
end1=cell2mat(end1);
end2=cell2mat(end2);
end1=sign(end1(1,:)).*min(abs(end1));
end2=sign(end2(1,:)).*max(abs(end2));
ends=[end1 end2];

function xx=maxndim(x,dims)
xx=x;
for i_dim=1:length(dims)
    xx=max(xx,[],dims(i_dim));
end
