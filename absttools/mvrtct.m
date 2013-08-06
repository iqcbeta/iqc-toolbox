function Z=mvrtct(X,Y)
% function Z=mvrtct(X,Y)
%
% forms a tf lti system Z=[Z;Y]
%
% Written by ameg@mit.edu to by-pass the bug in the
% Control Systems Toolbox v. 5.1
% Last modified by cmj on 2013/4/18

[ax,bx,cx,dx]=ssdata(ss(X));
[ay,by,cy,dy]=ssdata(ss(Y));

if size(dx,2)~=size(dy,2),
    disp_str(6)
end

nx=size(ax,1);
ny=size(ay,1);
mx=size(dx,1);
my=size(dy,1);

a=[ax zeros(nx,ny);zeros(ny,nx) ay];
b=[bx;by];
c=[cx zeros(mx,ny);zeros(my,nx) cy];
d=[dx;dy];

Z=ss(a,b,c,d);
