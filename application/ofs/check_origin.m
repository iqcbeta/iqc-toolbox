function pas=check_origin(ay,Q0,Q1,tmax2)
%   check_t2(Pi'*a*Pi,Q0,Q1,tmax2)  Checks if stability conditions
%   are met for t>tmax2
%   i.e.,  Q1-Ft2'*Q0*Ft2>0   for t>tmax2?    (For d=0)
%
%   Version: 0.1    Last modified: March 3, 2001
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%

%
% change call of "value" to "value_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu


n=size(ay,1)+1;
Q=care(ay,zeros(n-1),eye(n-1),eye(n-1),zeros(n-1),eye(n-1));

abst_init_lmi;
lmitbx_options([0 0 0 0 1]);
lam=symmetric;
lam>0;
Q*ay+ay'*Q+lam*Q<0;
lmi_mincx_tbx(-lam);
lam=value_iqc(lam);

abst_init_lmi;
lmitbx_options([0 0 0 0 1]);
k0=symmetric;
k0*Q-Q0>0;
lmi_mincx_tbx(k0);
k0=value_iqc(k0);

pas=0;
if min(eig(Q1-k0*exp(-lam*tmax2)*Q))>0
  pas=1;
end

