function [f,Y,Z]=iqc_white(n,b,a)
% function [f,Y,Z]=iqc_white(n,b,a)
%
% IQCs for deterministic white signals over the frequency intervall [-b,b].
%
% a contains the poles of the multiplier
% n denotes the number of components of the signal
%
% The IQCs are of the form
%
% f'*Y*f>0
%             N
% where Y=X0+sum Zk/(s+a(k))+conj(Zk)/(s+conj(a(k)))
%            k=1
%
% subject to the constraint
%             N
%    tr[b*X0+sum re[Zk*arctan(b/a(k))]]>0
%            k=1
% We use Zk=Xk+iYk, where Xk,Yk are real nxn matrices
%
% The output Z is a cell array containing the Xk and Yk in the order
%
%   [X1 Y1]
% Z=[  :  ]
%   [XN YN]
%   [X0  0]
%
% Inputs: n size of signal
%         b bandwidth (default=10)
%         a vector with pole locations (default a=1)
%
% Work by Ulf Jonsson October 1997
% Last modified by cmj on 2013/5/5

if nargin<1
    n=1;
end
if nargin<2
    b=10;
end
if nargin<3
    a=1;
end

global ABST

switch ABST.systemtype
    case 'continuous'
        s=tf([1 0],1);
        f=signal(n);
        Y=OO(n);
        Cstr=OO(1);
        N=length(a);
        for k=1:N
            if imag(a(k))==0
                Z{k,1}=rectangular(n); %#ok<*AGROW>
                Z{k,2}=OO(n);
                Y=Y+Z{k,1}*(1/(s+a(k)));
                Cstr=Cstr+atan(b/a(k))*trace(Z{k,1});
            else
                Z{k,1}=rectangular(n);
                Z{k,2}=rectangular(n);
                A=[-2*real(a(k))*eye(n) -abs(a(k))*eye(n);eye(n) zeros(n,n)];
                B=[eye(n);zeros(n,n)];
                C=[Z{k,1} real(a(k))*Z{k,1}+imag(a(k))*Z{k,2}];
                D=zeros(2*n,n);
                Ysl=ss(A,B,eye(2*n),D);
                Y=Y+C*Ysl;
                Cstr=Cstr+real(atan(b/a(k)))*trace(Z{k,1})-imag(atan(b/a(k)))*trace(Z{k,2});
            end
        end;
        Z{N+1,1}=rectangular(n);
        Z{N+1,2}=OO(n);
        Y=Y+Z{N+1,1};
        Cstr=Cstr+b*trace(Z{N+1,1});
        Cstr>0;
        f'*Y*f>0;
    case 'discrete'
        disp_str(70,'iqc_white','discrete')
end
