%% create velocity fields as MAT files with 0 on ground instead of as NC files

addpath(genpath('/home/orika/'))
disp('added path')

% monthly

%fid=fopen('soda_filenames_ALL.txt');
fid = fopen('soda_filenames_monthly_all.txt');
tline = fgetl(fid);
filename = cell(0,1);
while ischar(tline)
    z= strfind(tline , 'soda');
    filename{end+1,1} = tline(z(2):end) ;
    tline = fgetl(fid);
end
fclose(fid);

L = length(filename);
disp('created file list')

[U, V, W, longitude, latitude, depth] = load_fixed_velocity_fields_cluster(filename{1}, 1);
time = ncread(filename{1}, 'time');
savename = [filename{1}(1:end-3) '_MATfile.mat'];
save(savename, 'U', 'V', 'W', 'time');
save('griddata.mat', 'longitude', 'latitude', 'depth');
finished = 1

for i = 2:L
    filenamenow = filename{i};
    [U, V, W, ~, ~, ~] = load_fixed_velocity_fields_cluster(filenamenow, 1);
    time = ncread(filenamenow, 'time');
    savename = [filenamenow(1:end-3) '_MATfile.mat']
    save(savename, 'U', 'V', 'W', 'time');
    finished = finished + 1
end

