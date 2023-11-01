function [un] = EOM_polarMIT(x,u)
Lon = x(:,1);
Lat = x(:,2);
[LonPerM,LatPerM] = LonLatPerM(Lon,Lat);

U=u(:,1);
V=u(:,2);
W=u(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert values in ux and uy to dLon/dSec and dLat/dSec,
% respectively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = U.* LonPerM;
V = V.* LatPerM;

un(:,1)=U;
un(:,2)=V;
un(:,3)=W;
end