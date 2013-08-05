function sat(G,d,nu1,nu2a,nu2b);
%OFS  Global stability of equilibrium points of saturation systems.
%
%   SAT(G,D,NU1,NU2a,NU2b) proves global asymptotic stability
%   of the origin of system G in feedback with an on_off controller.
%   It proves the 3 impact maps associated with the system are
%   quadratically stable.  Using the notation in the paper, this 
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
%   By default are set to 30.
%
%   For example
%
%       s=tf([1 0],1);
%       G=-1.5*(2*s^2+2*s+1)*(s-1)/((s^3+2*s^2+6*s+2)*(s+5));
%	sat(G,1);
%
% For more information and related papers go to
%
% http://web.mit.edu/jmg/www/
%
%   Version: 1.0 		    Last modified: March 3, 2001
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%

ttt=which('d');
if isempty(ttt)
        d=1;
end
clear ttt

ttt=which('nu1');
if isempty(ttt)
        nu1=30;
end
clear ttt

ttt=which('nu2a');
if isempty(ttt)
        nu2a=nu1;
end
clear ttt

ttt=which('nu2b');
if isempty(ttt)
        nu2b=nu1;
end
clear ttt

ttt=which('npoi1');
if isempty(ttt)
        npoi1=300;
end
clear ttt

ttt=which('npoi2a');
if isempty(ttt)
        npoi2a=npoi1;
end
clear ttt

ttt=which('npoi2b');
if isempty(ttt)
        npoi2b=npoi1;
end
clear ttt


glob_sat(G,d,nu1,nu2a,nu2b,npoi1,npoi2a,npoi2b);
