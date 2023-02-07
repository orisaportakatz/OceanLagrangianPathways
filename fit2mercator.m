
function [U] = fit2mercator(u,lon)
% U=fit2mercator(u,lon);
% u- 4D velocity field from mitgcm
% first dim is longitude
%
% builed the variable u such as it will go from -180:180
% instead from 0:360. e.g grinich is in the middle of the map.
% given the already fited lon from -180 until 180 of u
%
% reason: to fit to a mecator projection

[~, in]=min(abs(lon-0));

U=cat(1,u(in+1:end,:,:,:),u(1:in,:,:,:));


