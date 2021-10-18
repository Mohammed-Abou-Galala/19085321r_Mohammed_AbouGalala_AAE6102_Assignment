%--------------------------------------------------------------------------
%   Single POINT Positioning (SPP) using  pseudorange AND phase observations 
%  ------------------------------------------------------------------------ 
%  Coder : Mohammed Abou-Galala
%  Date  : 13-10-2021 
%--------------------------------------------------------------------------
function [az,el,dist]=SatElAz(x1,y1,z1,phi,plam,x2,y2,z2)
%--------------------------------------------------------------------------
% syntax:
%   [az,el,dist]=SatElAz(x1,y1,z1,phi,plam,x2,y2,z2)
%
% input:        x1       = x-coordinate for the receiver position (m).
%               y1       = y-coordinate for the receiver position (m).
%               z1       = z-coordinate for the receiver position (m).
%               phi      = the latitude for the receiver position. 
%               lam      = the longitude for the receiver position. 
%               x2       = x-coordinate for the satellite position (m).
%               y2       = y-coordinate for the satellite position (m).
%               z2       = z-coordinate for the satellite position (m).
%
% output:       az       = azimuth 
%               el       = elevation angle 
%               dist     = the geometric distance 
%               idot     = rate of inclination angle (r/s)
%
% description : calculate the receiver-satellite geometric distance,
% elevation angle and azimuth. 
%--------------------------------------------------------------------------
rlat=phi-pi/2.0;
rlon=plam-pi;
srlat=sin(rlat);
crlat=cos(rlat);
srlon=sin(rlon);
crlon=cos(rlon);

dx=x2-x1;
dy=y2-y1;
dz=z2-z1;

du=crlat*crlon*dx+ crlat*srlon*dy -srlat*dz;
dv=srlon*dx      - crlon*dy;
dw=srlat*crlon*dx+ srlat*srlon*dy+ crlat*dz;

dist= sqrt(du^2+dv^2+dw^2);
az  = atan2(dv,du)/pi*180;
el  = asin(dw/dist)/pi*180;