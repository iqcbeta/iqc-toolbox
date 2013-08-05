function [w,p,q,b,c]=iqc_delay(v,T,a)
% function [w,p,q,b,c]=delay(v,T,a)
% 
% IQC description of the uncertain delay block
% w(t)=v(t-t0), 0<t0<T  (v can be a vector)
% 
% q*(|v|^2-|w|^2)>0 for ANY frequency-dependent weight q
%
% p|(p1p/q1p)(w+v')|^2>
%    p|(w-v')/q1m|^2-k1m*p|(w-v')/q1m|^2
%   +(2/T)*p*|(v-w)/q2m|^2-(2/T)*k2m*p|(v-w)/q2m|^2,
% where Re(p)>0,
%     |p1p/q1p|^2 is the upper bound of r1(W)
%     (1-k1m|s|^2)/|q1m|^2 is the lower bound of r1(W)
%     (1-k2m|s|^2)/|q2m|^2 is the lower bound of r2(W)
% W=wT/2,
%       r1(W)=sin(W)/W for W<pi
%            =0 otherwise
%       r2(W)=cos(W) for W<pi,
%            =0 otherwise
% q=b{1}+b{2}/(s+a(1))+b{2}/(s+a(2))+...
% p=c{1}+c{2}/(s+a(1))+c{3}/(s+a(2))+...
%
% v must be differentiable, though it is possible to remove
% this condition by re-writing the block
%
% Written by ameg@mit.edu, last modified Dec. 07, 1997

s=tf([1 0],1);      % for convenience, s=d/dt
m=size(v,1);        % input size
w=signal(m);        % define the output

% coefficients
T1=T/2;
T2=T1^2;
T3=T1^3;
T5=sqrt(2/T);
p1p=[T2*0.064 0 1];    % coefficients of filters
q1p=[T3*0.0292 T2*0.1957 T1*0.6553 1];
k1m=1/pi^2;
q1m=[T2*0.06037466827221   T1*0.43138708785584   1];
k2m=0.4073;
q2m=[T2*0.09219544457293   T1*0.52639423357960   1];

% dc term of q(s)
bb=rectangular(m); 
b{1}=bb;
q=bb;
v'*bb*v>w'*bb*w;

% basic signals
vd=derivative(v);
[y1p,x1p]=lti_full(p1p,q1p,w+vd);
[y2m,x2m]=lti_full(T5,q2m,v-w);
[y1m,x1m]=lti_full(1,q1m,w-vd);

cc=rectangular(m);
p=cc;
c{1}=cc;
y2md=derivative(y2m);
y1md=derivative(y1m);
% one IQC expressed using three inequalities
y1p'*cc*y1p>0;
y1m'*cc*y1m<y1md'*k1m*cc*y1md;
y2m'*cc*y2m<y2md'*k2m*cc*y2md;

for k=1:length(a),
   aa=a(k);
   Ga=abst(1/(s+aa));
   % define necessary signals
   wa=Ga*w;
   va=Ga*v;
   vad=derivative(va);
   ya1p=lti_a(p1p,q1p,aa,x1p,wa+vad);
   ya1m=lti_a(1,q1m,aa,x1m,wa-vad);
   ya2m=lti_a(T5,q2m,aa,x2m,va-wa);
   % q(s) terms
   bb=rectangular(m);
   b{k+1}=bb;
   q=q+bb*Ga;
   v'*bb*va>w'*bb*wa;
   % p(s) terms
   cc=rectangular(m);
   p=p+cc*Ga;
   c{k+1}=cc;
   % one IQC expressed using three inequalities
   ya1p'*cc*y1p>0;
   ya2md=derivative(ya2m);
   ya1md=derivative(ya1m);
   ya1m'*cc*y1m<ya1md'*k1m*cc*y1md;
   ya2m'*cc*y2m<ya2md'*k2m*cc*y2md;
end
p>0;
q>0;    