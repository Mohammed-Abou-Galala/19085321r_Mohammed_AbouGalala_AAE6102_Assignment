%--------------------------------------------------------------------------
%   Single POINT Positioning (SPP) using  pseudorange AND phase observations 
%  ------------------------------------------------------------------------ 
%  Coder : Mohammed Abou-Galala
%  Date  : 13-10-2021 
%--------------------------------------------------------------------------
% there are FOUR scenarios for SPP: 
% ---------------------------------
% options:
% observations     :  using only codes    or using codes and phases together.
% cutoff angle     :  using zero degree   or using 10 degree.
% residual limit   :  using 1000 m        or using 10 m.
% tropsphere       :  ignore tropo        or estimate the tropo.
% weights          :  using equal weights or using elevation weights.
%==========================================================================
%                            SCENARIO FOUR	
%==========================================================================
% the scenario four using the following options:
% observations     :  using  codes and phases.
% cutoff angle     :  using 10 degree.
% residual limit   :  using 10 m.
% tropsphere       :  estimate the troposphere.
% weights          :  using elevation weights.
% ambiguities      : estimate the ambiguities.
%==========================================================================
% The difference betwen the scenario one and scenario five is the limit for 
% the residual check is 10 meter rather than 1000 m in the scenario one,
% and the cutoff elevation angle is 10 degree rather than zero degree, and
% model the tropshere and estimate the wet troposphere,
% and using the elevation weights, and
% using the phase and codes together.
% we used here the single frequency ionosphere-free combination in order to
% elliminate the ionosphere effect.
%==========================================================================
clear;
%--------------------------------------------------------------------------
% initialized vectors and matrices
xsat=zeros(8,6); dtsat=zeros(8,1); dist=zeros(8,1); el=dist; az=dist; 
%--------------------------------------------------------------------------
% constants
clight = 299792458;         % Speed of light. (m/s)
Wedot= 7.2921151467e-5;     % WGS 84 value of earth’s rotation rate (r/s).
mu= 3.986005e+14;           % WGS 84 value of earth's universal gravitation constant (m^3/s^2)
F= -4.442807633e-10;        % Relativistic correction term constant.
lam1=clight/1575.42e6;
%--------------------------------------------------------------------------
% files
file_nav  ='eph.dat';
file_obs  ='rcvr.dat';
file1     = 'satposfile_scen4.txt'; fid_file1 = fopen(file1,'w');
file2     = 'rcvposfile_scen4.txt'; fid_file2 = fopen(file2,'w');
file3     = 'satresfile_scen4.txt'; fid_file3 = fopen(file3,'w');
fprintf(fid_file1,'PRN              X(m)                     Y(m)                    Z(m)                    SAT CLOCK(sec)  \n');
fprintf(fid_file2,'ITERATION         DX(m)        DY(m)          DZ(m)      DCLOCK(m)      DTROP(m)   FINAL X(m)      FINAL Y(m)       FINAL Z(m)    FINAL CLOCK(m) FINAL TROP(m)');
fprintf(fid_file2,'  SIG X(m)  SIG Y(m) SIG Z(m)  SIG CLOCK(m) SIG TROP(m) DOP\n');
fprintf(fid_file3,'ITERATION      PRN     AZIMUTH      ELEV     CODE RESIDUAL(m) \n');
m = 4;                          % number of unknowns or estimated parameters
%--------------------------------------------------------------------------
% read satellite navigation data
[tow_orb,iprn_orb,toe,a0,a1,a2,e,roota,dn,m0,idot,i0,omega,omega0,omegadot,cus,cuc,cis,cic,crs,crc] =read_broadcast_info(file_nav);
%--------------------------------------------------------------------------
% read observation data
[tow,iprn_obs,code,phase] =read_obs(file_obs);
n = length(code);   %number of pseudorange observations or equations
%--------------------------------------------------------------------------
% compute satellite positioins and clocks 
for i=1:length(iprn_obs)
    iprn = iprn_obs(i);
    index = find(iprn_orb==iprn);
    [xsat(i,1:3), dtsat(i,1)] = satposs((tow(index,1)),code(index,1),toe(index,1),...
        a0(index,1),a1(index,1),a2(index,1), toe(index,1), m0(index,1),roota(index,1)^2,...
        dn(index,1),e(index,1),omega(index,1), omega0(index,1), omegadot(index,1),i0(index,1),...
        idot(index,1),cus(index,1),cuc(index,1),cis(index,1),cic(index,1),crs(index,1),crc(index,1),toe(index,1));
end
%-------------------------------------------------------------------------
% print the satellite positions & clocks 
for i=1:length(iprn_obs)
    iprn = iprn_obs(i);
    fprintf(fid_file1,'G%02d       %+15.3f          %+15.3f          %+15.3f          %+20.15f\n',iprn, xsat(i,1),xsat(i,2),xsat(i,3),dtsat(i,1));   
end
fclose(fid_file1);
%--------------------------------------------------------------------------
% initial receiver coordinates/receiver clock offset
rece_pos   = [-2694685.473; -4293642.366; 3857878.924];
clock_bias = 0;
Cxx        = zeros(5,5);
Cxx(1,1)=1/(10^20); Cxx(2,2)=Cxx(1,1); Cxx(3,3)=Cxx(1,1); Cxx(4,4)=Cxx(1,1);
Cxx(5,5)=1/(0.01^2);
%--------------------------------------------------------------------------
% geodetic coordinates
[wlat, wlon, walt] = Wgsxyz2lla(rece_pos);
%--------------------------------------------------------------------------
% intial zenith troposphere delay
zen_trop   =  tropo_error_correction(90, walt);
%--------------------------------------------------------------------------
% initial ambiguity
for isat=1:length(iprn_obs)
    amb(isat) = (code(isat,1)-lam1*phase(isat,1))/lam1;
end
% loop for code observations
% ------------------------------------------------------------------------
iter=1;
while (1)
    while(1)
        icount=1; A=[]; computed_range=[]; observed_range=[]; Q=[];
        az=[];  el=[]; dist=[];
        for isat=1:length(iprn_obs)
%--------------------------------------------------------------------------
% elevation/azimuth/distace computation
            [az(isat,1),el(isat,1),dist(isat,1)]=SatElAz(rece_pos(1),rece_pos(2),rece_pos(3),wlat*pi/180,wlon*pi/180,xsat(isat,1),xsat(isat,2),xsat(isat,3));
%------------------------------------------------------------------------
% design or coefficient matrix
            if (el(isat,1)>10)
                A(icount,1) = (rece_pos(1) - xsat(isat,1))/dist(isat,1); % coefficient for X coordinate unknown
                A(icount,2) = (rece_pos(2) - xsat(isat,2))/dist(isat,1); % coefficient for Y coordinate unknown
                A(icount,3) = (rece_pos(3) - xsat(isat,3))/dist(isat,1); % coefficient for Z coordinate unknown
                A(icount,4) =  1;                                        % coefficient for receiver clock OFFSET 
                A(icount,5) =  1/sind(el(isat,1));                       % coefficient for tropo
                A(icount,5+icount) =  +0.5*lam1;                         % coefficient for ambiguity
%--------------------------------------------------------------------------
% troposphere modelling        
                corr = tropo_error_correction(el(isat,1), walt);
% -------------------------------------------------------------------------
% computed_range vector
                computed_range(icount,1) = dist(isat,1) + clock_bias - clight*dtsat(isat,1) + 1/sind(el(isat,1))* zen_trop - 0.5*lam1*amb(isat); 
% ------------------------------------------------------------------------
% observed range vector
                observed_range(icount,1) = 0.5*code(isat,1)+0.5*lam1*phase(isat,1) ;
% ------------------------------------------------------------------------
% observation covariance matrix
                  Q(icount,icount) = 1 /(sin(el(isat,1) * pi/180)^2);    % using the elevation weights.  
                  Cxx(5+icount,5+icount) = 1 /(2^2);
%                 Q(icount,icount) = 1 ;        
                  icount=icount+1;            
            end            
        end
% ------------------------------------------------------------------------        
        % remove data for low elevation
        for k=length(iprn_obs):-1:1
            if (el(k)<10)
                [iprn_obs,xsat,dtsat,code,phase,amb]=reject_obs_phase(iprn_obs,xsat,dtsat,code,phase,amb,iprn_obs(k));                
            end            
        end
% ------------------------------------------------------------------------        
        % misclosure vector
        diff = observed_range-computed_range;
% ------------------------------------------------------------------------
% least squares adjustment 
        Weight=diag((diag(Q).^-1));
% ------------------------------------------------------------------------
        % normal matrix
        N = (A'*(Weight)*A)+Cxx;
        if (length(iprn_obs)>= m) 
            % least squares solution
            x   = (N^-1)*A'*(Weight)*(observed_range-computed_range);
            % computation of residuals 
            residuals = (A*x + computed_range)-observed_range;        
        end
        
        % check residuals
        [maxres,ireject,iprn_rej] = findmax_residuals(residuals,iprn_obs);
        if (maxres > 10)
            [iprn_obs,xsat,dtsat,code,phase,amb,Cxx]=reject_obs_phase_res(iprn_obs,xsat,dtsat,code,phase,amb,Cxx,iprn_rej);
        else
            break;
        end
    end
% ------------------------------------------------------------------------    
    % update the solution
    rece_pos   = rece_pos   + x(1:3);
    clock_bias = clock_bias + x(4);
    zen_trop   = zen_trop   + x(5);
    
    for z=1:length(amb)
        amb(z) = amb(z)+x(5+z);        
    end
% ------------------------------------------------------------------------    
    % check the convergence 
    dx = norm(x(1:3));   
% ------------------------------------------------------------------------    
    % a posteriori covariance matrix of the estimation error
    posteriori_sigma = (residuals'*(Weight)*residuals)/(icount-m);
    Cxx = (N^-1);
    cov_XR  = Cxx(1:3,1:3);
    var_dtR = Cxx(4,4);
    std_x  = sqrt(cov_XR(1,1));  std_y  = sqrt(cov_XR(2,2));   std_z  = sqrt(cov_XR(3,3));
    std_dtR = sqrt(var_dtR);    std_trop = sqrt(Cxx(5,5));
% ------------------------------------------------------------------------    
    % DOP
    cov_XYZ = (A'*A)^-1;
    PDOP   = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
% ------------------------------------------------------------------------   
    % print pos results
    fprintf(fid_file2,'%04d         %+12.5f  %+12.5f  %+12.5f  %+11.3f  %+12.5f',iter,x(1),x(2),x(3),x(4),x(5));   
    fprintf(fid_file2,'  %+12.3f    %+12.3f    %+12.3f    %+12.3f %+12.3f',rece_pos(1),rece_pos(2),rece_pos(3),clock_bias,zen_trop);  
    fprintf(fid_file2,'    %+8.3f  %+8.3f  %+8.3f  %+8.3f    %+8.3f     %+5.2f\n',std_x,std_y,std_z,std_dtR,std_trop, PDOP);
% ------------------------------------------------------------------------        
    % print residuals 
    for j=1:length(iprn_obs)
        iprn_p = iprn_obs(j); index=find(iprn_obs==iprn_p);
        if (az(index)<0), az(index)=az(index)+180; end
        fprintf(fid_file3,'%02d              G%02d     %6.2f      %5.2f         %+6.2f  \n',iter,iprn_p,az(index),el(index),residuals(index));
    end
    fprintf(fid_file3,'                                           \n');
% ------------------------------------------------------------------------
    % check the convergence 
    if (dx<1e-4), break; end
    iter=iter+1;
    if (iter>10), break; end         % the maximum number of iterations 10. 
end
fclose(fid_file2); fclose(fid_file3);
fclose('all');