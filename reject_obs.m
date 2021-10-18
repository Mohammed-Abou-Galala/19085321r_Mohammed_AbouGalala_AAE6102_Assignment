%--------------------------------------------------------------------------
%   Single POINT Positioning (SPP) using  pseudorange AND phase observations 
%  ------------------------------------------------------------------------ 
%  Coder : Mohammed Abou-Galala
%  Date  : 13-10-2021 
%--------------------------------------------------------------------------
function [iprn_obs,xsat,dtsat,code]=reject_obs(iprn_obs,xsat,dtsat,code,iprnrej)
%--------------------------------------------------------------------------
% syntax:
%   [iprn_obs,xsat,dtsat,code]=reject_obs(iprn_obs,xsat,dtsat,code,iprnrej)
%
% input:        iprn_obs = vector for the observed satellite.
%               xsat     = matrix for satellite positions. 
%               dtsat    = vector for satellite clock corrections. 
%               code     = vector for the observed codes.
%               iprnrej  = the iprn of the rejected satellite. 
%
% output:       iprn_obs = vector for the observed satellite.
%               xsat     = matrix for satellite positions. 
%               dtsat    = vector for satellite clock corrections. 
%               code     = vector for the observed codes.
%
% description : remove the rejected code residual. 
%--------------------------------------------------------------------------

index = find(iprn_obs==iprnrej);

iprn_obs(index)=[];
xsat(index,:)=[];
dtsat(index)=[];
code(index)=[];
% phase(index)=[];