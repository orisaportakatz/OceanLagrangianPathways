function fullMonthly_reg_run_cluster(startmonth, startyear, Nyears, savename, configname, initialconditions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demo_3D_GlobalOcean_SODA_MOWtraj.m
%
% Particle tracking in the North Atlantic, specifically from the MOW region
% 3D Velocity field from SODA
%
% Yael Amitai ; Israela Musan 27.2.2019 clean ; Ori S Katz 3.3.2019 ;
% Nadav Mantel 29.11.2019
% Ori Saporta-Katz 9/05/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input variables:

%startmonth/startyear = date of start of pathways calculation
%Nyears = number of years to calculate pathways
%savename = filename to save the trajectories
%configname = name of config file
%initialconditions = matrix of initial conditions, syntax:
%initialconditions(longitude, latitude, depth) such that
%size(intialconditions) = (number-of-initial-conditions, 3).

%% configuration

addpath(genpath('/home/orika/'));

run(configname); %configure parameters.

% system("ps -u orika")  %debug

%startup;

%% file list
%create a list of all MAT file we want to work with, soda_filename.txt contain
% all data since 1980 to 2018s:
fid=fopen('soda_filenames_monthly_all.txt');
tline = fgetl(fid);
filename_nc = cell(0,1);
while ischar(tline)
    z= strfind(tline , 'soda');
    filename_nc{end+1,1} = tline(z(2):end) ;
    tline = fgetl(fid);
end
fclose(fid);

L = length(filename_nc);

for i = 1:L
    filename{i} = [filename_nc{i}(1:end-3) '_MATfile.mat'];
end

startyear_index = startyear - 1980 + 1;
runUntil = startyear_index + Nyears - 1;

%% initialization:

%% load grid and velocity data, fixed to include only North Atlantic with W regridded to fit U and V's domain:
load('griddataMonthly.mat', 'longitude', 'latitude', 'depth'); %already cut to include only North Atlantic
filename_1 = filename{startyear_index};
%[U, V, W, longitude, latitude, depth] = load_fixed_velocity_fields_cluster(filename_1, 1);
load(filename_1, 'U', 'V', 'W', 'time');
U1 = U;
V1 = V;
W1 = W;
time1 = time; %time in days since 1/1/1980 at 0:00:00
disp('loaded first velocity field')

%if issteady or isperiodic == 1, the velocity field needs to be created only once.
%In these cases, the next lines create 3 13-month velocity fields: U, V, W

if issteady %when issteady==1, the run is on a steady flow repeating startmonth/startyear over and over again. Note that this induces particles getting stuck due to small eddies close to borders!
    U = repmat(U(:,:,:,startmonth),1,1,1,13);
    V = repmat(V(:,:,:,startmonth),1,1,1,13);
    W = repmat(W(:,:,:,startmonth),1,1,1,13);
    regrun = 0;
    tsnapshot = 1;
    disp('created relevant velocity field')
    
elseif isperiodic %when isperiodic==1, repeat startyear over and over again, starting from startmonth.
    U = cat(4, U(:,:,:,startmonth:end), U(:,:,:,1:startmonth));
    V = cat(4, V(:,:,:,startmonth:end), V(:,:,:,1:startmonth));
    W = cat(4, W(:,:,:,startmonth:end), W(:,:,:,1:startmonth));
    regrun = 0;
    disp('created relevant velocity field')
else %if issteady==0 and isperiodic==0 this is a regular time-dependent run
    regrun = 1;
end

%% specify boundary conditions:
if regrun == 0
    if isfreeslip
        [U2, V2] = UVfreeSlip(U, V, 3); %inserts free slip condition and turns all ground points to nan!
        U = U2;
        V = V2;
    end
    if ispushShore
        [U3, V3] = UVpushShore(U, V, 3);
        U = U3;
        V = V3;
    end
    %make sure ground is 0 and not nan:
    U(isnan(U)) = 0;
    V(isnan(V)) = 0;
    W(isnan(W)) = 0;
    
    %% for 2d flow:
    if is2dim
        W = zeros(size(W));
    end
end

%% grid:
[Lon,Lat,Depth]=ndgrid(longitude,latitude,depth);
gridxyz=cell(1,3);% see notbook.pdf for explanation
gridxyz{1}=Lon;
gridxyz{2}=Lat;
gridxyz{3}=Depth;

%% Initial conditions
% Create initial conditions:
if ~exist('initialconditions', 'var')
    m_9 = 200; %number of intial points in latitude
    n_9 = 200; %number of intial points in depth
    
    YIC = linspace(34, 37, m_9);
    %ZIC = linspace(0, 1000, n_9);
    ZIC = 618.7031;
    XIC = -7;
    [XICmesh, YICmesh, ZICmesh] = meshgrid(XIC, YIC, ZIC);
    X = XICmesh(:);
    Y = YICmesh(:);
    Z = ZICmesh(:);
    
    cc = zeros(length(X), 3);
    cc(:,1) = X; cc(:,2) = Y; cc(:,3) = Z;
else
    cc = initialconditions;
end

xyzAll = reshape(cc,1,length(cc(:,1)),3);

% system("ps -u orika")  %debug

%% iteration over all years
for y = startyear_index:runUntil
    
    %% create velocity fields for regular run:
    %if this is a regular run (regrun==1, issteady==0 and isperiodic==0),
    %the next lines create 3 13-month velocity fields: U, V, W - fitting
    %years y through y+1:
    if regrun == 1
        filename_1 = filename{y}
        filename_2 = filename{y+1}
        %create a 13-month year from startmonth/currentyear to startmonth/currentyear+1:
        if y == startyear_index
            U1 = U(:,:,:,startmonth:end);
            V1 = V(:,:,:,startmonth:end);
            W1 = W(:,:,:,startmonth:end);
            time1 = time(startmonth:end);
        else
            U1 = U_next_year(:,:,:,startmonth:end);
            V1 = V_next_year(:,:,:,startmonth:end);
            W1 = W_next_year(:,:,:,startmonth:end);
            time1 = time(startmonth:end);
        end
        %[U_next_year, V_next_year, W_next_year, ~, ~, ~] = load_fixed_velocity_fields_cluster(filename_2, 1);
        load(filename_2, 'U', 'V', 'W', 'time');
        U_next_year = U;
        V_next_year = V;
        W_next_year = W;
        time_next_year = time;
        
        U = cat(4,U1,U_next_year(:,:,:,1:startmonth));
        V = cat(4,V1,V_next_year(:,:,:,1:startmonth));
        W = cat(4,W1,W_next_year(:,:,:,1:startmonth));
        %debug: whos U
        
        tsnapshot = [time1;time_next_year(1:startmonth)]; % time vector fitting the velocity vector snapshots, signifying time in days for U,V,W since 1/1/1980 at 0:00:00.
        tsnapshot_sec = tsnapshot * 24 * 60 * 60;
        %debug: tsnapshot
        
        % system("ps -u orika")  %debug
        
        %% specify boundary conditions:
        isfreeslip
        if isfreeslip
            [U2, V2] = UVfreeSlip(U, V, 3); %inserts free slip condition and turns all ground points to nan!
            U = U2;
            V = V2;
        end
        if ispushShore
            [U3, V3] = UVpushShore(U, V, 3);
            U = U3;
            V = V3;
        end
        
        %make sure ground is 0 and not nan:
        U(isnan(U)) = 0;
        V(isnan(V)) = 0;
        W(isnan(W)) = 0;
        
        %% for 2d flow:
        if is2dim
            W = zeros(size(W));
        end
        disp(['created relevant velocity field - year ' num2str(1980 + y - 1)])
        
        % system("ps -u orika")  %debug
        
    end
    
    %% Initialization:
    if y == startyear_index
        
        %Create a bathymetry mask for the velocity - 0 in ground and 1 in ocean:
        GroundIndicator = sign(U(:,:,:,1).^2 + V(:,:,:,1).^2 + W(:,:,:,1).^2).^2;
        
        if ismask
            %interpolated mask - 1 in ocean bulk, 0 in ground, linear
            %interpolation between ocean and ground on "shore" area
            Umask_prep = griddedInterpolant(Lon, Lat, Depth, sign(GroundIndicator), 'linear', 'none');
            if issharpmask % if ==1, make mask sharp - 1 in ocean bulk and shore, 0 in ground
                Umask = @(x,y,z) 1 - (Umask_prep(x,y,z) == 0);
            else
                Umask = Umask_prep;
            end
        else %if no mask, use this as a degenerate mask, just for syntax
            Umask = @(x,y,z) 1;
        end
        % system("ps -u orika")  %debug
        
        
        %new initial conditions for the next year:
    elseif (y == runUntil) %at the last step:
        cc = [squeeze(xyz(end,:,1))' squeeze(xyz(end,:,2))' squeeze(xyz(end,:,3))'];
    else %set initial conditions as last position of trajectories as calculated for previous iteration
        cc = [squeeze(xyz(end,:,1))' squeeze(xyz(end,:,2))' squeeze(xyz(end,:,3))'];
        % system("ps -u orika")  %debug
    end
    
    %% time settings:
    
    timetoolbox = (tsnapshot_sec - tsnapshot_sec(1)); %time in sec, timetoolbox(1)=0
    %whos timetoolbox %debug
    
    %Here we dump every 2.5 days because SODA defines first 10 months as 30
    %days and last 2 months as 32.5 days
    tdump=24*60*60 *2.5;  % (1 day in sec * x days)
    
    %timeindex = 1;
    isfirstime = 1;
    direction = 1; %forward
    %         direction=-1; %backward
    
    y %output which year index we are in
    
    %% yearly loop:
    disp('starting yearly loop')

    for t=1:12 %for the given year "y", loop over each month+succesive month, starting from startmonth
        %tic  %count how long a month takes
        
        % Prepare for tracking
        u=cell(1,3);% see notbook.pdf for explanation
        u{1}=U(:,:,:,t:t+1);
        u{2}=V(:,:,:,t:t+1);
        u{3}=W(:,:,:,t:t+1).*-1; % positive W is directed up;
        
        tspan = timetoolbox(t):tdump:timetoolbox(t+1); % save position at these time steps
        tspan_offset = tsnapshot_sec(t):tdump:tsnapshot_sec(t+1);
        
        %Here I switch config to use ode15s to solve bug3 maybe?? --Ori 8.7.19
        %             [tmp_xyz,tmp_ts] = lpta_ode_outer(gridxyz, u, cc,'config_GlobalOcean_SODA.m',TT(t:t+1),tspan);
        %[tmp_xyz,tmp_ts] = lpta_ode_outer_new(gridxyz, u, cc,configname,TT(t:t+1),tspan,Umask);
        [tmp_xyz,tmp_ts] = lpta_ode_outer_new(gridxyz, u, cc,configname,timetoolbox(t:t+1),tspan,Umask);
        
        % system("ps -u orika")  %debug
        
        %whos tmp_ts
        %whos tmp_xyz
        cc = [ tmp_xyz(end,:,1)', tmp_xyz(end,:,2)' , tmp_xyz(end,:,3)'];
        
        %% insert integration output into entire array:
        if (direction==1)
            
            if isfirstime
                xyz = tmp_xyz;
                
                isfirstime = 0;
                %timevec = tmp_ts + tsnapshot_sec(1);
                
                %whos tmp_xyz
                
            else
                xyz = cat(1, xyz, tmp_xyz(2:end, : , :));
                %timevec = [timevec; tmp_ts(2:end) + timevec(end) - timevec(1)];
                %whos tmp_xyz
            end
        end
        
        clear u
       
    end
    
    disp(['finished yearly loop for year ', num2str(startyear + y - 1)])
    
    xyzAll = cat(1,xyzAll, xyz(2:end,:,:));
    % system("ps -u orika")  %debug
end

xyzAll = xyzAll(1:10:end,:,:); %undersampling result to 25 day time-steps for memory conservation.
whos xyzAll
save(savename, 'xyzAll');

disp('the end')
