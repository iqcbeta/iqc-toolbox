%%%%%% Collection of systems
%
%   Version: 1.1    Last modified: Aug 9, 1999
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%

s=tf([1 0],1)

%  3rd order

G=tf([-1 -1 4],[1 6 11 6]) 	% original one
G=tf([-4 -4 4],[1 4 6 4])	%   "	    "  with complex poles
G=tf([1],[1 3 3 1])		% relative order 3   1/(s+1)^3
G=(s^2+3*s+10)/((s^2+4*s+2)*(s+3))      % CB=1
G=(s^2-0.02*s+0.98)/((s^2+0.02*s+1.02)*(s+.1))


%  5h order

G=0.5*tf([-3 -0.5 7.5 12.75 7], [1 4.5 8.75 8.75 4.625 1.125])  % orginal
G=-(s+1)^2*(s-3)*(s+2)/((s^2+2*s+3)*(s+2)*(s+10)^2)

%  6th order

G=-(s+1)^2*(s-3)*(s+2)^2/((s^2+2*s+3)*(s+2)*(s+10)^2*(s+1))
G=(s+1)*(s^2+4*s+10)*(s-5)^2/((s^2+2*s+2)*(s^2+5*s+5)^2)  % CB=1 and y goes up
G=(s+1)*(s^2-4*s-10)*(s+4*s-5)/((s^2+2*s+2)*(s^2+5*s+12)^2) % CB=0, CAB=5. g has 3 zeros, but only one is a lc. Need Cones.


%  7th order

G=-(s+1)^2*(s-3)*(s+2)^3/((s^2+2*s+3)*(s+2)*(s+10)^2*(s^2+10*s+1))
G=-(s^2+2*s+1)*(s-3)*(s+2)^3/((s^2+2*s+3)*(s+2)*(s+10)^2*(s^2+10*s+1))


%%%%%%%%%

