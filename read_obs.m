%--------------------------------------------------------------------------
%   Single POINT Positioning (SPP) using  pseudorange AND phase observations 
%  ------------------------------------------------------------------------ 
%  Coder : Mohammed Abou-Galala
%  Date  : 13-10-2021 
%--------------------------------------------------------------------------
function [tow,iprn,code,phase] =read_obs(file)
%--------------------------------------------------------------------------
% syntax:
%   [tow,iprn,code,phase] =read_obs(file)
%
% input:        file =  file for observation data.
%
% output:       tow      = time of week(s)
%               iprn     = satellite PRN number
%               code     = code observations (m).
%               phase    = carrier phase observations (m). 

% description : read the code and phase observations. 
%--------------------------------------------------------------------------
% initialize        
tow =zeros(8,1);
iprn=zeros(8,1);
code =zeros(8,1);
phase  =zeros(8,1);
frac   =zeros(8,1);

isat=1;

fid=fopen(file,'rt');

while(~feof(fid))
    
    line=fgetl(fid);
    if (length(line)>1)
        str=sscanf(line(1:16),'%c');
        tow(isat,1)=str2double(str);
        str=sscanf(line(18:19),'%c');
        iprn(isat,1)=str2double(str);
        str=sscanf(line(25:36),'%c');
        code(isat,1)=str2double(str);
        str=sscanf(line(43:49),'%c');
        phase(isat,1)=str2double(str);
        str=sscanf(line(51:54),'%c');
        frac(isat,1)=str2double(str);
        phase(isat,1)=phase(isat,1)+(frac(isat,1)/2048);
        
        isat=isat+1;    
    end
end
fclose(fid);