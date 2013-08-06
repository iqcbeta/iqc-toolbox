function Z=mhrzct(X,Y)
% function Z=mhrzct(X,Y)
%
% forms a tf lti system Z=[Z Y]
%
% Written by ameg@mit.edu to by-pass the bug in the
% Control Systems Toolbox v. 5.1
% Last modified by cmj on 2013/4/18

[ax,bx,cx,dx]=ssdata(ss(X));
[ay,by,cy,dy]=ssdata(ss(Y));

if size(dx,1)~=size(dy,1),
    disp_str(6)
end

nx=size(ax,2);
ny=size(ay,2);
mx=size(dx,2);
my=size(dy,2);

a=[ax zeros(nx,ny);zeros(ny,nx) ay];
b=[bx zeros(nx,my);zeros(ny,mx) by];
c=[cx cy];
d=[dx dy];

Z=ss(a,b,c,d);
