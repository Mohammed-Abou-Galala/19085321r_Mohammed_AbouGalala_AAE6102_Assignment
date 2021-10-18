%--------------------------------------------------------------------------
%   Single POINT Positioning (SPP) using  pseudorange AND phase observations 
%  ------------------------------------------------------------------------ 
%  Coder : Mohammed Abou-Galala
%  Date  : 13-10-2021 
%--------------------------------------------------------------------------
function [tow,iprn,toe,a0,a1,a2,e,roota,dn,m0,idot,i0,omega,omega0,omegadot,cus,cuc,cis,cic,crs,crc] =read_broadcast_info(file)
%--------------------------------------------------------------------------
% syntax:
%   [tow,iprn,toe,a0,a1,a2,e,roota,dn,m0,idot,i0,omega,omega0,omegadot,cus,cuc,cis,cic,crs,crc] =read_broadcast_info(file)
%
% input:        file =  file for navigation data.
%
% output:       tow      = time of week(s)
%               iprn     = satellite PRN number
%               toe      = reference time of ephemeris parameters (s)
%               a0       = clock correction coefficient – group delay (s)
%               a1       = clock correction coefficient (s/s)
%               a2       = clock correction coefficient (s/s/s)
%               e        = eccentricity
%               roota    = square root of semi-major axis a (m**1/2)
%               dn       = mean motion correction (r/s)
%               m0       = mean anomaly at reference time (r)
%               idot     = rate of inclination angle (r/s)
%               i0       = inclination angle at reference time (r)
%               omega    = argument of perigee (r)
%               omegadot = right ascension (r)
%               cus      = argument of latitude correction, sine (r)
%               cuc      = argument of latitude correction, cosine (r)
%               cis      = inclination correction, sine (r)
%               cic      = inclination correction, cosine (r)
%               crs      = radius correction, sine (m)
%               crc      = radius correction, cosine (m)

% description : read the information for the satellite positions and
% satellite clock correction.
%--------------------------------------------------------------------------
% initialize        
tow =zeros(8,1);
iprn=zeros(8,1);
toe =zeros(8,1);
a0  =zeros(8,1);
a1  =zeros(8,1);
a2  =zeros(8,1);
e   =zeros(8,1);
roota =zeros(8,1);
dn    =zeros(8,1);
m0    =zeros(8,1);
idot  =zeros(8,1);
i0    =zeros(8,1);
omega    =zeros(8,1);
omega0   =zeros(8,1);
omegadot =zeros(8,1);

cus =zeros(8,1);
cuc =zeros(8,1);
cis =zeros(8,1);
cic =zeros(8,1);
crs =zeros(8,1);
crc =zeros(8,1);

isat=1;

fid=fopen(file,'rt');

while(~feof(fid))
    
    line=fgetl(fid);
    if (length(line)>1)
        str=sscanf(line(1:16),'%c');
        tow(isat,1)=str2double(str);
        str=sscanf(line(18:19),'%c');
        iprn(isat,1)=str2double(str);
        str=sscanf(line(41:53),'%c');
        toe(isat,1)=str2double(str);
        str=sscanf(line(57:70),'%c');
        a0(isat,1)=str2double(str);
        str=sscanf(line(74:87),'%c');
        a1(isat,1)=str2double(str);
        str=sscanf(line(92:104),'%c');
        a2(isat,1)=str2double(str);
        str=sscanf(line(134:146),'%c');
        e(isat,1)=str2double(str);
        str=sscanf(line(159:171),'%c');
        roota(isat,1)=str2double(str);
        str=sscanf(line(184:196),'%c');
        dn(isat,1)=str2double(str);
        str=sscanf(line(208:221),'%c');
        m0(isat,1)=str2double(str);
        str=sscanf(line(233:246),'%c');
        omega(isat,1)=str2double(str);
        str=sscanf(line(258:271),'%c');
        omega0(isat,1)=str2double(str);
        str=sscanf(line(284:296),'%c');
        i0(isat,1)=str2double(str);
        str=sscanf(line(300:313),'%c');
        omegadot(isat,1)=str2double(str);
        str=sscanf(line(317:330),'%c');
        idot(isat,1)=str2double(str);
        
        str=sscanf(line(335:347),'%c');
        cus(isat,1)=str2double(str);
        str=sscanf(line(351:364),'%c');
        cuc(isat,1)=str2double(str);
        str=sscanf(line(368:381),'%c');
        cis(isat,1)=str2double(str);
        str=sscanf(line(385:398),'%c');
        cic(isat,1)=str2double(str);
        str=sscanf(line(402:415),'%c');
        crs(isat,1)=str2double(str);
        str=sscanf(line(420:432),'%c');
        crc(isat,1)=str2double(str);
        
        isat=isat+1;    
    end
end
fclose(fid);