function [U, V, W, longitude, latitude, depth] = load_fixed_velocity_fields_cluster(filename, isNC, vel_variable)
%input:
%"filename" is a filename from which longitude, latitude, and depth are derived.
%If isNC = 1,  "filename" is used to derive also U, V, W, and no need to
%specify "vel_variable".
%If isNC = 0, "vel_variable" is used to derive U, V, W. Usually used for climatological annual and monthly. 

if isNC == 1
    U = ncread(string(filename),'u');
    V = ncread(string(filename),'v');
    W = ncread(string(filename),'wt');
else
    U = vel_variable.u;
    V = vel_variable.v;
    W = vel_variable.wt;
end
    
%load longitude, latitudes, depth:
longitude = ncread(string(filename),'xu_ocean'); % 1 (0degE) to 720 (360degE)
latitude = ncread(string(filename),'yu_ocean'); % 1 (-90degN) to 330 (90degN)
depth = ncread(string(filename),'sw_ocean'); % 1 (+5.0335 m) to 50 (+5395 m) - positive z direction is _downward_, so positive w velocity is _downwelling_.

longitude = longitude - 180; %set from -180 to 180 for mercator projection

%Prepare U:

U(isnan(U)) = 0;
U = fit2mercator(U,longitude);

%continue working on U:
U = U(1:end-1, 1:end-1, :, :);
% Reduce vel field to North Atlantic:
U = U(120:410, 129:310,:,:);

%Now, prepare V:

V(isnan(V)) = 0;
V = fit2mercator(V,longitude);
V = V(1:end-1, 1:end-1, :, :);
% Reduce vel field to North Atlantic:
V = V(120:410, 129:310,:,:);

%Finally, prepare W:

W(isnan(W)) = 0;
W = fit2mercator(W,longitude);

% %Regridding W to same grid as U, V: original W set for xt_ocean, yt_ocean
% %instead of xu_ocean, yu_ocean.
W = 1/2 * (W(:,1:end-1,:,:) + W(:, 2:end, :, :));
W = 1/2 * (W(1:end-1,:,:,:) + W(2:end, :, :, :));

% Reduce vel field to North Atlantic:
W = W(120:410, 129:310,:,:);

%now, turn into zeros all point in W for which U and V are ground:
W = W .* (sign(U).^2);
W(isnan(W)) =  0;

% Next, shorten longitude and latitude accordingly:
longitude = longitude(1:end-1);
latitude = latitude(1:end-1);
% %Reduce scope to North Atlantic
latitude = latitude(129:310);
longitude = longitude(120:410);








