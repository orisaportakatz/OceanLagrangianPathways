%bug3 - Ori - 16.7.19 - patch 1/7 - add "signOcean" to inner input parameters.
function [xyz,ts] = lpta_ode_inner_new(gridxyz, u, t2pts,tspan, xyz0, interpolant, EOM ,K,dt_stoch,options, odesolver, signOcean , Umask)

%interpolant = 'linear';

% syntax : [xyz,ts] = lpta_ode_inner(gridxyz, u, t2pts,tspan, xyz0, interpolant, EOM ,K,dt_stoch,options, odesolver)
% PARTICLE_TRACK_ODE_INNER  Generates particle tracks from a set of currents
% defined on a general grid using ODE solver
%
% This function is a (very) modified version of David M. Kaplan, which was inspired by Bruce Lipphardt's trajectories code.
% the main difference is the fact that we keep the interpolation coefficients, instead of re-calculating them in each time step of
% the ode solver.
% Inputs
% ------
% ---- required inputs:
% - gridxyz = is a cell array of dimensions (u-components,x-components) ;
% if the grid is non-staggered, it could be (1, x-components)
% each cell in gridxyz contains a matrix with dimensions (x,y,z), denoting the positions of the velocity samples;
% 3 dimensions is just an example. one could have 1 or 2 as well.
% - u         = is a cell array of dimensions (1,u-components), of the gridded velocity samples. each cell is of size (nx,ny,nz,2)
% - t2pts     = time points corresponding to the the given velocity samples
% - tspan     = times of dumps. in the case of a steady field run
% - xyz0      = tracer's initial conditions; it is a matrix (not cell array) with dimensions (n_tracers,n_dimensions)
% - interpolant - a string, typically 'spline' or 'linear', to instruct our algorithm
%           how to interpolate between grid points. If the u fields include NaN's, 'spline' can't be used.
% - EOM       = a function handle, representing the equation of motion $dx/dt=EOM(x,u)$.
%                      In most cases, @EOM(x,u)=u. see documentation for further info
% - K         = momentum eddy diffusivity. Either a scalar, or a vector with number of elements = number of dimensions.
% - dt_stoch  = time step for the Euler Maruyama scheme
% - options   = a data structure produced by matlab's [[http://www.mathworks.com/help/matlab/ref/odeset.html][odeset]].
% - odesolver = a function handle to your own ode solver, or one of matlab's (e.g. @ode45)
% - signOcean = predefined sign of particles in the ocean. Default is +1,
% i.e. positive Z means the particle is inside the ocean. This is defined
% in the new config files.
% - Uind is created with velocity fields, it is of the same size as the
%velocity field for a single time-step, with 0 in ocean and 1 in ground.


%
% OUTPUTS:
%
% xyz = 3 dimensional matrix with dimensions [time tracer_index components], indicating the positions of the tracer in the dumping times;
% ts  = list of dump times
%----------------------------------------------------------------------------------------------------
% Author : 2015 Avi Gozolchiani
% $Log: lpta_ode_inner.m,v $
% Revision 1.2  2015/03/09 15:58:33  avigoz
% adding time range (tbounds) and copyright headers
%
% Revision 1.7  2015/03/09 15:47:18  avigoz
% Summary: including tbounds - beginning and end time of the particles
%
%----------------------------------------------------------------------------------------------------

%% defaults
if ~exist( 'odesolver', 'var' )
    odesolver = @ode45;
end

if ~exist( 'options', 'var' )
    options = [];
end
nsteps=2;
if(length(t2pts)~=nsteps)
    error('t2pts should include 2 elements');
end
if(t2pts(1)==t2pts(2)) % steady state case
    nsteps=1;
end


% we don't check the input since it's have already been checked in the
% outer iteration
n_components=length(u);

%% prepare interpolants
% calculate the interpolation coefficients, before time stepping
for i_component=1:n_components % by component we mean velocity component, not space vectors
    for i_stp=1:nsteps
        if(n_components==1)
            tmp=griddedInterpolant(gridxyz{i_component,1},u{i_component}(:,i_stp),interpolant);
            %uinterp{i_component,i_stp}=@(xyz)tmp(xyz(:,1));
            %uinterp{i_component,i_stp} =@(xyz)tmp(xyz(:,1),xyz(:,2),xyz(:,3)) .* Umask(xyz(:,1),xyz(:,2),xyz(:,3));
            uinterp{i_component,i_stp} =@(xyz)tmp(xyz(:,1),xyz(:,2),xyz(:,3));
        elseif(n_components==2)
            tmp=griddedInterpolant(gridxyz{i_component,1},gridxyz{i_component,2},u{i_component}(:,:,i_stp),interpolant);
            %uinterp{i_component,i_stp}=@(xyz)tmp(xyz(:,1),xyz(:,2));
            %uinterp{i_component,i_stp} =@(xyz)tmp(xyz(:,1),xyz(:,2),xyz(:,3)) .* Umask(xyz(:,1),xyz(:,2),xyz(:,3));
            uinterp{i_component,i_stp} =@(xyz)tmp(xyz(:,1),xyz(:,2),xyz(:,3));
        elseif(n_components==3)
            
            %%%%%%%%%%%
            %bug3 - Ori - 16.7.19 - patch 2/7 - add "nearest" at the end of griddedInterpolant.
            %%% This defines the extrapolator as taking the value of the velocity
            %%% field of the nearest neighbor instead of extrapolating with the interpolant option.
            tmp=griddedInterpolant(gridxyz{i_component,1},gridxyz{i_component,2},gridxyz{i_component,3},u{i_component}(:,:,:,i_stp),interpolant,'nearest');
            %Old version:
            %tmp=griddedInterpolant(gridxyz{i_component,1},gridxyz{i_component,2},gridxyz{i_component,3},u{i_component}(:,:,:,i_stp),interpolant);
            %%%%%%%%%%%
            
            %     uinterp{i_component,i_stp}=@(xyz)tmp(xyz(:,1),xyz(:,2),xyz(:,3));
            uinterp{i_component,i_stp} =@(xyz)tmp(xyz(:,1),xyz(:,2),xyz(:,3)) .* Umask(xyz(:,1),xyz(:,2),xyz(:,3));
            %uinterp{i_component,i_stp} =@(xyz)tmp(xyz(:,1),xyz(:,2),xyz(:,3));
                        
        else % should never get here
            error('we currently don''t support more than 3D');
        end
    end % for i_stp=1:nsteps
end % for i_component=1:n_components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do tracking.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntraj = size(xyz0,1);
y0=xyz0(:);

if(nsteps==1)
    tfrom=tspan(1);
    tto=tspan(end);
else
    tfrom=t2pts(1);
    tto=t2pts(2);
end
if(K(1)<=0)
    dt_stoch=tto-tfrom;
end % if(K<=0)
ts=[];
xyz=[];
while(tfrom<tto)
    cur_t2pts=[tfrom tfrom+dt_stoch];
    cur_tspan=tspan(tspan>=cur_t2pts(1) & tspan<=cur_t2pts(2));
    cur_tspan2=unique([cur_t2pts(1) cur_tspan(:)' cur_t2pts(2)]); % merge dump points and the end points
    
    %%%%%%%%%%%
    %bug3 - Ori - 16.7.19 - patch 3/7: add signOcean to arguments of odesolver
    %change 6:
    %[tmpts,A] = odesolver(@ptrack_ode_worker, cur_tspan2, y0, options, ...
    %    uinterp,cur_t2pts,n_components,EOM ,nsteps==1, signOcean);
    longitude = squeeze(gridxyz{1,1}(:,1,1));
    latitude = squeeze(gridxyz{1,2}(1,:,1));
    depth = squeeze(gridxyz{1,3}(1,1,:));
    
    longitude = longitude(:);
    latitude = latitude(:);
    depth = depth(:);
    tic
    [tmpts,A] = odesolver(@ptrack_ode_worker, cur_tspan2, y0, options, ...
        uinterp,cur_t2pts,n_components,EOM ,nsteps==1, signOcean);
    toc
    %%%%%%%%%%%
    
    %%%%%%%%%%%
    %bug3 - Ori - 16.7.19 - patch 4/7 - implement the vertical boundary condition on A (instead of on tmpxyz):
    if n_components==3 %vertical boundary condition - this sets the z component of the particles to 0 if particles leave the basin to the air, i.e. obtain a sign different from signOcean.
        Atemp = reshape(A, length(tmpts), ntraj, n_components);
        Az = Atemp(:,:,3);
        Az(find(sign(Az) ~= signOcean)) = 0;
        Atemp(:,:,3) = Az;
        A = reshape(Atemp, length(tmpts), ntraj*n_components);
    end
    %%%%%%%%%%%
    
    tfrom=tmpts(end);
    
    
    if(K(1)>0)
        y0=reshape(A(end,:),ntraj,n_components);
        y0=y0+repmat(sqrt(2*K(:)'*dt_stoch),ntraj,1).*randn(size(y0));
        y0=y0(:);
        
        %%%%%%%%%%%
        %bug3 - Ori - 16.7.19 - patch 5/7 - in case the odesolver crashes and the
        %loop starts again, call this out as an error and break:
    else
        if tfrom < tto
            error('The integration did not complete. Check particles at the boundaries.')
        end
        %%%%%%%%%%%
    end
    
    [~,subset_t,~]=intersect(tmpts,cur_tspan);
    
    tmpxyz=reshape(A(subset_t,:),[],ntraj,n_components);
    
    %bug3 - Ori - remove this:
    % if n_components==3 %vertical boundry condition - this sets the z component of the particles to 0 if particles leave the basin to the air, i.e. obtain a positive z value.
    %       tmpD=tmpxyz(:,:,3);
    %       lastDsignIn=find(sign(xyz0(:,3)),1);
    %       tmpD(tmpD~=0 & sign(tmpD)~=sign(xyz0(lastDsignIn,3)))=0;
    %       tmpxyz(:,:,3)=tmpD;
    %   end
    %
    
    
    tmpts=tmpts(subset_t);
    % concatenate the results of different snapshots
    if(~isempty(tmpts))
        xyz=cat(1,xyz,tmpxyz);
        ts=cat(1,ts,tmpts);
    end
end
[ts,its,~]=unique(ts);
xyz=xyz(its,:,:);
%%%%%%----------------------Subfunctions--------------------%%%%%%%
%% core

%%%%%%%%%%%
%bug3 - Ori - 16.7.19 - patch 6/7: add signOcean to arguments of ptrack
%change 5:
%function f = ptrack_ode_worker( T, y, uinterp,t2pts,n_components,EOM,flag_steady, signOcean)
function f = ptrack_ode_worker( T, y, uinterp,t2pts,n_components,EOM,flag_steady, signOcean)
%%%%%%%%%%%

%function f = ptrack_ode_worker( T, y, uinterp,t2pts,n_components,EOM,flag_steady,signAir)
% Break up X, Y and Z coordinates.
ntraj = length(y)/n_components;
xyz = reshape(y,ntraj,n_components);

%%%%%%%%%%%
%bug3 - Ori - 16.7.19 - patch 7/7: take particles that left to the air back to the surface
if n_components == 3
    ind = find(sign(xyz(:,3)) ~= signOcean);
    xyz(ind,3) = 0;
end
%%%%%%%%%%%


% see if we are in the relevant time bin
if(~flag_steady) % otherwise it's a steady field
    if(T-t2pts(1)<0 | t2pts(2)-T<0)
        f = nan(size(y));
        error('we are not in the relevant time bin');
    end % if(T-t2pts(1)<0 | t2pts(2)-T<0)
    % linearly interpolate in time
    % u(t)=(u1-u0)/dt*(t-t0)+u0 =  u1/dt*(t-t0)+u0*(1-(t-t0)/dt) = a1*u1+a0*u0
    % where a1=(t-t0)/dt , a0=(1-(t-t0)/dt)
    a1=(T-t2pts(1))/diff(t2pts);
    a0=1-a1;
end % if(t2pts(1)~=t2pts(2)) ... else ...

% make a mask for nans
nn = all(~isnan(xyz),2);
% interpolate u between 2 subsequent time snapshots
u=nan(ntraj,n_components);
for i_component=1:n_components
    if(flag_steady)
        %change 3:
        %u(nn,i_component)=uinterp{i_component,1}(xyz(nn,:));
        %u(nn,i_component)=uinterp{i_component,1}(xyz(nn,:)) .* landmask(xyz(:,1),xyz(:,2),xyz(:,3),Uind,longitude,latitude,depth);
        u(nn,i_component)=uinterp{i_component,1}(xyz(nn,:));
    else % if(flag_steady)
        %change 4:
        %u(nn,i_component)=(a0*uinterp{i_component,1}(xyz(nn,:))+a1*uinterp{i_component,2}(xyz(nn,:)));
        %u(nn,i_component)=(a0*uinterp{i_component,1}(xyz(nn,:))+a1*uinterp{i_component,2}(xyz(nn,:))) .* landmask(xyz(:,1),xyz(:,2),xyz(:,3),Uind,longitude,latitude,depth);
        u(nn,i_component)=(a0*uinterp{i_component,1}(xyz(nn,:))+a1*uinterp{i_component,2}(xyz(nn,:)));
    end % if(flag_steady) .. else ..
end
% our equation of motion is just (f=) dx/dt=u
u(nn,:)=EOM(xyz(nn,:),u(nn,:));

%%%%%%%%%%%
%bug3 - Ori - 16.7.19 - patch 8/8: take velocity of surface particles to zero.
%if n_components == 3
%    u(ind,3) = 0;
%end
%%%%%%%%%%%


f = u(:);
