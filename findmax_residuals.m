%--------------------------------------------------------------------------
%   Single POINT Positioning (SPP) using  pseudorange AND phase observations 
%  ------------------------------------------------------------------------ 
%  Coder : Mohammed Abou-Galala
%  Date  : 13-10-2021 
%--------------------------------------------------------------------------
function [maxres,ireject,iprn_rej] = findmax_residuals(res,iprn_obs)
%--------------------------------------------------------------------------
% syntax:
%   [maxres,ireject,iprn_rej] = findmax_residuals(res,iprn_obs)
%
% input:        res      = code residuals (m).
%               iprn_obs = vector for the observed satellite.
%
% output:       maxres   = maximun code residual. 
%               ireject  = index for the rejected residual.  
%               iprn_rej = the prn of the rejected satellite.  
%
% description : find the maximum code residual. 
%--------------------------------------------------------------------------
maxres=0; ireject=[];

for i=1:length(res)
    if (abs(res(i))>maxres)
       maxres   = abs(res(i));
       ireject  = i;
       iprn_rej = iprn_obs(i);
    end   
end