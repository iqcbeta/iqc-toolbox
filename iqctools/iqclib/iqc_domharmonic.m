function [f,M,x,d]=iqc_domharmonic(n,a,b,N1,N2)
% function [f,M,x,d]=iqc_domharmonic(n,a,b,N1,N2)
%
% IQCs for signals with dominant harmonics, i.e., the
% spectrum is concentrated to be in a certain frequency band.
% There are three alternatives, 1) bandpass, 2) lowpass, and
% 3) highpass characteristic:
%
% 1) supp f(jw)=[-b,-a]U[a,b]  (supp means support)
%
% 2) supp f(jw)=[-a,a], if b=[], a>0
%
% 3) supp f(jw)=(-inf,-|a|]U[|a|,inf), if b=[], a<0
%
% The IQCs are of the form
%
%  f'*(x*H'*H-d)*f>0,
%
%  where
%
%  1) H is N1 order butterworth BP filter
%     with cut off frequencies at w1=a/N2 and w2=b*N2, and
%     x*real(H(ja))-d>0
%     x*real(H(jb))-d>0
%     x,d>0
%
%  2) H is N1 order butterworth LP filter
%     with cut off frequency at w=a*N2, and
%     x*real(H(ja))-d>0
%     x,d>0
%
%  3) H is N1 order butterworth LP filter
%     with cut off frequency at w=|a|/N2, and
%     x*real(H(j|a|))-d>0
%     x,d<0
%
%   Inputs: n size of f, default n=1
%           a lower frequency bound, default a=1
%           b upper frequency bound, default b=[]
%           N1 order of butterworth filter, default N1=1
%           N2 determines break frequency of filter, default N2=2
%
%  Work by Ulf Jonsson, June 1998
%
if nargin<5, N2=2; end
if nargin<4, N1=1; end
if nargin<3, b=[]; end
if nargin<2, a=1; end
if nargin<1, n=1; end

if exist('butter')==0
    error('Your computer doesn''t have Signal Processing Toolbox!')
end

f=signal(n);
s=tf([1 0],1);
x=symmetric;
d=symmetric;
if isempty(b)
    if a>0   %case 2
        [num,den]=butter(N1,N2*a,'s');
        H=tf(num,den);
        x>0;
        d>0;
        x*real(freqresp(H'*H,a))-d>0;
    else     %case 3
        a=abs(a);
        [num,den]=butter(N1,a/N2,'s');
        H=tf(num,den);
        x<0;
        d<0;
        x*real(freqresp(H'*H,a))-d>0;
    end
else       %case 1
    a=abs(a);
    [num,den]=butter(N1,[a/N2,b*N2],'s');
    H=tf(num,den);
    x>0;
    d>0;
    x*real(freqresp(H'*H,a))-d>0;
    x*real(freqresp(H'*H,b))-d>0;
end
(H*f)'*x*(H*f)>f'*d*f;
M=H'*x*H-d;





