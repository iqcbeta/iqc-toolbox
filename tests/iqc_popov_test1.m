function iqc_popov_test1(n,cs)
% test for iqc_popov on the system
% v=G(s)(f-w),  w=d(v), d within sector [0,1]
% estimating the gain f->w
%
% the result is compared with the lower bound
% which is the maximum (over frequency) of the function
% 
% a(w)=|G/Im(G)| when Re(G)<-1
%     =|G/(G+1)| otherwise
%
% G of order n is generated randomly (default G=k(1-s)/(s^2+2s+5),
% cs is the case variable (default cs=1)
%
%
% Written by ameg@mit.edu,  last modified October 13, 1997

% default values
if nargin<1,
   G=tf([-1 1],[1 2 5]);
else
   if nargin<2,
      cs=1;
   end
   randn('state',cs);
   A=randn(n);  A=A-(1+max(real(eig(A))))*eye(n);
   B=randn(n,1);
   C=randn(1,n);
   G=ss(A,B,C,0);
end
G=G/(norm(G,inf)+1);

% calculating the lower bound for the gain
[m,p,ww]=bode(G);
ww=linspace(0,max(ww)/3,1000);
g=freqresp(G,ww);
g=squeeze(g);
lbw=abs(g./(1+g)).*(real(g)>-1)+abs(g./imag(g)).*(real(g)<=-1);
lb=max(lbw);


abst_init_iqc
f=signal;
w=signal;
v=G*(f-w);
w==iqc_popov(v,[0 1]);
iqc_gain_tbx(f,w)
disp(['The lower bound is ' num2str(lb)])
iqc_bode
