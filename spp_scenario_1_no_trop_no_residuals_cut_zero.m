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
%                            SCENARIO ONE
%==========================================================================
% the scenario one using the following options:
% observations     :  using only codes.
% cutoff angle     :  using zero degree.
% residual limit   :  using 1000 m.
% tropsphere       :  ignore troposphere and ionosphere.
% weights          :  using equal weights.
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
%--------------------------------------------------------------------------
% files
file_nav  ='eph.dat';
file_obs  ='rcvr.dat';
file1     = 'satposfile_scen1.txt'; fid_file1 = fopen(file1,'w');
file2     = 'rcvposfile_scen1.txt'; fid_file2 = fopen(file2,'w');
file3     = 'satresfile_scen1.txt'; fid_file3 = fopen(file3,'w');
fprintf(fid_file1,'PRN              X(m)                     Y(m)                    Z(m)                    SAT CLOCK(sec)  \n');
fprintf(fid_file2,'ITERATION     DX(m)        DY(m)          DZ(m)      DCLOCK(m)     FINAL X(m)      FINAL Y(m)     FINAL Z(m)    FINAL CLOCK(m)');
fprintf(fid_file2,'   SIG X(m)  SIG Y(m)  SIG Z(m)  SIG CLOCK(m) DOP\n');
fprintf(fid_file3,'ITERATION      PRN     AZIMUTH      ELEV     CODE RESIDUAL(m) \n');
m = 4;                       % number of unknowns or estimated parameters
%--------------------------------------------------------------------------
% read satellite navigation data
[tow_orb,iprn_orb,toe,a0,a1,a2,e,roota,dn,m0,idot,i0,omega,omega0,omegadot,cus,cuc,cis,cic,crs,crc] =read_broadcast_info(file_nav);
%--------------------------------------------------------------------------
% read observation data
[tow,iprn_obs,code,phase] =read_obs(file_obs);
n = length(code);            % number of pseudorange observations or equations
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
% Cxx        = zeros(4,4);
%--------------------------------------------------------------------------
% geodetic coordinates
[wlat, wlon, walt] = Wgsxyz2lla(rece_pos);
%--------------------------------------------------------------------------
% loop for code observations
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
            if (el(isat,1)>0)
                A(icount,1) = (rece_pos(1) - xsat(isat,1))/dist(isat,1); % coefficient for X coordinate unknown
                A(icount,2) = (rece_pos(2) - xsat(isat,2))/dist(isat,1); % coefficient for Y coordinate unknown
                A(icount,3) = (rece_pos(3) - xsat(isat,3))/dist(isat,1); % coefficient for Z coordinate unknown
                A(icount,4) =  1;                                        % coefficient for receiver clock OFFSET 
                %   A(icount,5) =  1/sind(el(isat,1));                   % coefficient for tropo
%--------------------------------------------------------------------------
% troposphere modelling        
                corr = tropo_error_correction(el(isat,1), walt);
% -------------------------------------------------------------------------
% computed_range vector
                computed_range(icount,1) = dist(isat,1) + clock_bias - clight*dtsat(isat,1); 
% ------------------------------------------------------------------------
% observed range vector
                observed_range(icount,1) = code(isat,1) ;
% ------------------------------------------------------------------------
% observation covariance matrix
%         Q(icount,icount) = 1 /(sin(el(isat,1) * pi/180)^2);
                Q(icount,icount) = 1 ;    %  equal weights       
                icount=icount+1;     
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
        N = (A'*(Weight)*A);
        if (length(iprn_obs)>= m) 
            % least squares solution
            x   = (N^-1)*A'*(Weight)*(observed_range-computed_range);
            % computation of residuals 
            residuals = (A*x + computed_range)-observed_range;        
        end
        
        % check residuals
        [maxres,ireject,iprn_rej] = findmax_residuals(residuals,iprn_obs);
        if (maxres > 1000)
            [iprn_obs,xsat,dtsat,code]=reject_obs(iprn_obs,xsat,dtsat,code,iprn_rej);
        else
            break;
        end
    end
% ------------------------------------------------------------------------    
    % update the solution
    rece_pos   = rece_pos   + x(1:3);
    clock_bias = clock_bias + x(4);
    % corr_zen   = corr_zen + x(5);
% ------------------------------------------------------------------------        
    % check the convergence 
    dx = norm(x(1:3));   
% ------------------------------------------------------------------------        
    % a posteriori covariance matrix of the estimation error
    posteriori_sigma = (residuals'*(Weight)*residuals)/(icount-m);
    Cxx = posteriori_sigma*(N^-1);
    cov_XR  = Cxx(1:3,1:3);
    var_dtR = Cxx(4,4);
    std_x  = sqrt(cov_XR(1,1));  std_y  = sqrt(cov_XR(2,2));   std_z  = sqrt(cov_XR(3,3));
    std_dtR = sqrt(var_dtR);
% ------------------------------------------------------------------------        
    % DOP
    cov_XYZ = (A'*A)^-1;
    PDOP   = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
% ------------------------------------------------------------------------        
    % print pos results
    fprintf(fid_file2,'%d        %+12.5f  %+12.5f  %+12.5f  %+11.3f   ',iter,x(1),x(2),x(3),x(4));   
    fprintf(fid_file2,'%+12.3f    %+12.3f    %+12.3f    %+12.3f',rece_pos(1),rece_pos(2),rece_pos(3),clock_bias);  
    fprintf(fid_file2,'    %+8.3f  %+8.3f  %+8.3f  %+8.3f     %+5.2f\n',std_x,std_y,std_z,std_dtR, PDOP);
% ------------------------------------------------------------------------        
    % print residuals 
    for j=1:length(iprn_obs)
        iprn_p = iprn_obs(j); index=find(iprn_obs==iprn_p);
        if (az(index)<0), az(index)=az(index)+180; end
        fprintf(fid_file3,'%d              G%02d     %6.2f      %5.2f         %+6.2f  \n',iter,iprn_p,az(index),el(index),residuals(index));
    end
    fprintf(fid_file3,'                                           \n');
% ------------------------------------------------------------------------        
    % check the convergence 
    if (dx<1e-4), break; end
    iter=iter+1;
    if (iter>10), break; end  % the maximum number of iterations 10. 
end
fclose(fid_file2); fclose(fid_file3);
fclose('all');