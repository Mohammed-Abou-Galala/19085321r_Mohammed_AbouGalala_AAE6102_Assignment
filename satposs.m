%--------------------------------------------------------------------------
%   Single POINT Positioning (SPP) using  pseudorange AND phase observations 
%  ------------------------------------------------------------------------ 
%  Coder : Mohammed Abou-Galala
%  Date  : 13-10-2021 
%  source: RTKLIB OPEN SOURCE 2007-2020 by T.TAKASU 
%--------------------------------------------------------------------------
function [rs, dts] = satposs(time,pr,toc,f0,f1,f2, toe, M0,A,deln,e,omg,OMG0,OMGd,i0,idot,cus,cuc,cis,cic,crs,crc,toes)
%--------------------------------------------------------------------------
% syntax:
%   [rs, dts] = satposs(time,pr,toc,f0,f1,f2, toe, M0,A,deln,e,omg,OMG0,OMGd,i0,idot,cus,cuc,cis,cic,crs,crc,toes)
%
% input:        time     = time of week(s)
%               pr       = satellite PRN number
%               toc      = reference time of clock parameters (s)
%               f0       = clock correction coefficient – group delay (s)
%               f1       = clock correction coefficient (s/s)
%               f2       = clock correction coefficient (s/s/s)
%               toe      = reference time of ephemeris parameters (s)
%               e        = eccentricity
%               A        = semi-major axis a (m)
%               deln     = mean motion correction (r/s)
%               M0       = mean anomaly at reference time (r)
%               idot     = rate of inclination angle (r/s)
%               i0       = inclination angle at reference time (r)
%               omg      = argument of perigee (r)
%               OMGd     = right ascension (r)
%               cus      = argument of latitude correction, sine (r)
%               cuc      = argument of latitude correction, cosine (r)
%               cis      = inclination correction, sine (r)
%               cic      = inclination correction, cosine (r)
%               crs      = radius correction, sine (m)
%               crc      = radius correction, cosine (m)
%
% output:       rs      = satellite posituion (m)
%               dts     = satellite clock correction includes relativeity correction

% description : calculate the satellite position and satellite clock coreection 
% using the broadcast ehemeris. 
%--------------------------------------------------------------------------
CLIGHT = 299792458.d0;
% transmission time by satellite clock 
time=time-(pr/CLIGHT);
% satellite clock bias by broadcast ephemeris 
dt = eph2clk(time, toc, f0, f1, f2);
time=time-dt;
% satellite position and clock at transmission time 
[rs, dts] = eph2pos(time,toe,M0,A,deln,e,omg,OMG0,OMGd,i0,idot,cus,cuc,cis,cic,crs,crc,toes,f0,f1,f2);
end 
% -------------------------------------------------------------------------
function eph2clk = eph2clk(time, toc, f0, f1, f2)
% -------------------------------------------------------------------------
%  input:       time     = time of week (sec) 
%               toc      = reference time of clock parameters (s)
%               f0       = clock correction coefficient – group delay (s)
%               f1       = clock correction coefficient (s/s)
%               f2       = clock correction coefficient (s/s/s)
%  output:      eph2clk  = satellite clock correction without relativeity correction
%  description : compute satellite clock bias with broadcast ephemeris (gps)
%  and satellite clock bias (s) without relativeity correction
% ----------------------------------------------------------------------------
t=time-toc;
for i=1:2
    t=t-(f0+f1*t+f2*t*t);
end
eph2clk=f0+f1*t+f2*t*t;
end
%--------------------------------------------------------------------------
function [rs, dts] = eph2pos(time,toe,M0,A,deln,e,omg,OMG0,OMGd,i0,idot,cus,cuc,cis,cic,crs,crc,toes,f0,f1,f2)
% -------------------------------------------------------------------------
% input:        time     = time of week(s)
%               toc      = reference time of clock parameters (s)
%               f0       = clock correction coefficient – group delay (s)
%               f1       = clock correction coefficient (s/s)
%               f2       = clock correction coefficient (s/s/s)
%               toe      = reference time of ephemeris parameters (s)
%               e        = eccentricity
%               A        = semi-major axis a (m)
%               deln     = mean motion correction (r/s)
%               M0       = mean anomaly at reference time (r)
%               idot     = rate of inclination angle (r/s)
%               i0       = inclination angle at reference time (r)
%               omg      = argument of perigee (r)
%               OMGd     = right ascension (r)
%               cus      = argument of latitude correction, sine (r)
%               cuc      = argument of latitude correction, cosine (r)
%               cis      = inclination correction, sine (r)
%               cic      = inclination correction, cosine (r)
%               crs      = radius correction, sine (m)
%               crc      = radius correction, cosine (m)
% description : compute satellite position and clock bias with broadcast ephemeris (gps)
% satellite clock includes relativity correction.
% -----------------------------------------------------------------------------
CLIGHT = 299792458.d0;
MU_GPS = 3.9860050d14;          % gravitational constant
OMGE   = 7.292115d-5;           % earth angular velocity (rad/s) ref (2) 
RTOL_KEPLER = 1d-13;            % relative tolerance for Kepler equation 
MAX_ITER_KEPLER = 30;           % max number of iteration of Kelpler 

if(A<=0.d0)
    rs=0.d0; dts=0.d0;
    return;
end

tk=(time-toe);
mu=MU_GPS; omge1=OMGE;
M=M0+(sqrt(mu/(A*A*A))+deln)*tk;
n=0;E=M;Ek=0.d0;

while(abs(E-Ek)>RTOL_KEPLER && n<MAX_ITER_KEPLER)
    Ek=E; E=E-(E-e*sin(E)-M)/(1.d0-e*cos(E));
    n=n+1;
end

if(n>=MAX_ITER_KEPLER), return; end
sinE=sin(E); cosE=cos(E);

u=atan2(sqrt(1.d0-e*e)*sinE,cosE-e)+omg;
r=A*(1.d0-e*cosE);
i=i0+idot*tk;
sin2u=sin(2.d0*u); cos2u=cos(2.d0*u);
u=u+cus*sin2u+cuc*cos2u;
r=r+crs*sin2u+crc*cos2u;
i=i+cis*sin2u+cic*cos2u;
x=r*cos(u); y=r*sin(u); cosi=cos(i);


O=OMG0+(OMGd-omge1)*tk-omge1*toes;
sinO=sin(O); cosO=cos(O);
rs(1)=x*cosO-y*cosi*sinO;
rs(2)=x*sinO+y*cosi*cosO;
rs(3)=y*sin(i);

tk=(time-toe);
dts=f0+f1*tk+f2*tk*tk;

% relativity correction 
dts=dts-2.d0*sqrt(mu*A)*e*sinE/(CLIGHT*CLIGHT);
end 
%--------------------------------------------------------------------------
