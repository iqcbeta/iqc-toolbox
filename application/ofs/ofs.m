function ofs(G,d,nu1,nu2,tmax1,tmax2,sproc);
%OFS  Global stability of equilibrium points of on_off systems.
%
%   OFS(G,D,NU1,NU2,TMAX1,TMAX2,SPROC) proves global asymptotic stability
%   of the equilibrium point of the system G in feedback with an on_off
%   controller.  It proves the 2 impact maps associated with the system
%   are quadratically stable.  Using the notation in the paper, this 
%   function finds the parameters of the quadratic Lyapunov functions
%   defined in the switching surface.  A graphic with the minimum
%   eigenvalues of
%               R_i(t)>0   for all ti-<= t_i <=ti+
%   for i=1,2, is plotted, where ti- and ti+ are, respectively,
%   bounds on the minimum and maximum expected switching times
%
%   G is an LTI system SYS (created with either TF or SS).
%
%   D is the parameter of the on_off controller.  By default is set to 1.
%
%   NUi are the numbers of LMIs used to sample the interval ti-<= t_i <=ti+.
%   By default are set to 40.
%
%   TMAXi are the maximum switching times where the LMIs will be checked.
%   Later verstions will find them automatically.  By default are set to 4.
%
%   SPROC if equal to 1 runs the less conservative condition based on the
%   S-Procedure.  By default is set to 0.
%
%   For example
%
%       s=tf([1 0],1);
%       G=-(2*s^2+2*s+12)/(s^3+2*s^2+2*s+3);
%	ofs(G,1);
%
% For more information and related papers go to
%
% http://web.mit.edu/jmg/www/
%
%   Version: 0.1	NOTE: this is an alpha version.
%   Last modified: March 3, 2001
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%

ttt=which('d');
if isempty(ttt)
        d=1;
end
clear ttt

ttt=which('nu1');
if isempty(ttt)
        nu1=40;
end
clear ttt

ttt=which('nu2');
if isempty(ttt)
        nu2=nu1;
end
clear ttt

ttt=which('tmax1');
if isempty(ttt)
        tmax1=4;
end
clear ttt

ttt=which('tmax2');
if isempty(ttt)
        tmax2=tmax1;
end
clear ttt

ttt=which('npoi1');
if isempty(ttt)
        npoi1=500;
end
clear ttt

ttt=which('npoi2');
if isempty(ttt)
        npoi2=npoi1;
end
clear ttt

ttt=which('sproc');
if isempty(ttt)
        sproc=0;
end
clear ttt

if d==0
  origin(G,nu1,nu2,tmax1,tmax2,npoi1,npoi2);
elseif sproc==1
  not_origin_c(G,d,nu1,nu2,tmax1,tmax2,npoi1,npoi2);
else 
  not_origin(G,d,nu1,nu2,tmax1,tmax2,npoi1,npoi2);
end


