function [w,Mrl]=iqc_ratelimiter(v,a,k)
% function [w,Mrl]=iqc_ratelimiter(v,a,k)
%
% Give IQCs for the "rate limiter" operator:
%
%         ________       _____       ________
%        |        | w1  |     |  x  |        |
% v ---->| sat(.) |---->| 1/s |---->| s+a(1) |--> w
%    -|  |________|     |_____|  |  |________|
%     |                          |
%     |__________________________|
%
%  where sat(.) is a semiconcave function with the
%  property d(sat(0))/dt=k.
%
%  IQCs given in this function is 
%  1. Sasha IQC for x, v-x and w1 
%  2. Zame-Falb IQCs of v-x and w1. Multipliers for ZF IQCs:
%     Xi/(s+a(i)), i=2,3,4...  
%
%  default value: a=[1,1],k=1
%
%  written by cykao@mit.edu  last modified: Feb 19. 1998
 
if nargin<3, k=1; end
if nargin<2, a=[1,1]; end
if nargin<1, error('input must be specified'), end
if size(v,1)>1, error('input must be scalar'), end

s=tf([1 0],1);
w=signal;
x=(1/(s+a(1)))*w;
w1=w-a(1)*x;

% Sasha's IQC for x, v-x, and w1
Mrl=symmetric(3);
m1=Mrl(1,1);  %% [ m1 m2 m3 ]
m2=Mrl(1,2);  %% | m2 m4 m5 |
m3=Mrl(1,3);  %% [ m3 m5 m6 ]
m4=Mrl(2,2);
m5=Mrl(2,3); 
m6=Mrl(3,3);
[m1,m2;m2,2*m4]>0;
m4+2*k*m5+k^2*m6>0;
m6<0;
[x;v-x;w1]'*Mrl*[x;v-x;w1]>0;


% Zames-Falb and Popov IQC for v-x and w1
vn=(500/(s+500))*(v-x);
w1==iqc_monotonic(vn,a(2:length(a)),0,k);

