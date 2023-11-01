function [LonPerM, LatPerM] = LonLatPerM( Lon, Lat, varargin )
% LONLATPERCM  Calculates the Longitude and Latitude displacement
% equivalent to 1 m at given locations on the earth.
%
% Usage: [LonPerM, LatPerM] = LonLatPerM( Lon, Lat, spheroid )
%
% This function uses m_map toolbox to make calculations
% 
% Inputs
% ------
% Lon,Lat = equal size matrices of Lon, Lat coordinates
% spheroid = see m_fdist for more details
%
% Outputs
% -------
% LonPerM = change in Lon for 1 cm movement in longitudinal direction at
%            locations in Lon,Lat.
% LatPerM = change in Lat for 1 cm movement in latitudinal direction at
%            locations in Lon,Lat.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: LonLatPerM.m 396 2009-01-15 16:39:00 efredj $	
%
% Copyright (C) 2009 Erick Fredj
% Licence: GPL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get XXX/M - remember angles in True convention
LonPerM = mod( m_fdist( Lon, Lat, 90, 1.0, varargin{:} ) - Lon, 360 );
[tt,ll] = m_fdist( Lon, Lat,  0, 1.0, varargin{:} );
LatPerM = mod( ll - Lat, 360 );
