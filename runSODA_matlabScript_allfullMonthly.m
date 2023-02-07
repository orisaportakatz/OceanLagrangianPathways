function runSODA_matlabScript_allfullMonthly(savenameprefix, configname, Nyears, runindex)

%ps = parallel.Settings;
%ps.Pool.AutoCreate = false;

%maxNumCompThreads(1);

addpath(genpath('/home/orika/'));

%% create year and month data:
allyears = 1980:1999;
allmonths = 1:12;
[AllYears, AllMonths] = meshgrid(allyears, allmonths);
AllYears = AllYears(:);
AllMonths = AllMonths(:);

startyear = AllYears(runindex);
startmonth = AllMonths(runindex);

%% create initial conditions for testing:
m_9 = 200; %number of intial points in latitude
n_9 = 201; %number of intial points in depth

YIC = linspace(33, 38, m_9);
ZIC = linspace(0, 1600, n_9);
%ZIC = 618.7031;
XIC = -7;
[XICmesh, YICmesh, ZICmesh] = meshgrid(XIC, YIC, ZIC);
X = XICmesh(:);
Y = YICmesh(:);
Z = ZICmesh(:);

cc = zeros(length(X), 3);
cc(:,1) = X; cc(:,2) = Y; cc(:,3) = Z;


    %% create savename
    savename = [savenameprefix '_fullMonthly_' num2str(startmonth, '%02.f') '_' ...
        num2str(startyear, '%02.f') '_' num2str(Nyears, '%02.f') 'yrs.mat'];
    %% run
startnow = 1
	fullMonthly_reg_run_cluster(startmonth, startyear, Nyears, savename, configname, cc)
end


