function [w,p,q,x,y,z]=iqc_ltigain(v,a)
% function [w,p,q,x,y,z]=iqc_ltigain(v,a)
%
% defines iqc's for the relation
% w(t)=delta*v(t),
% where delta is a constant scalar
% between -1 and 1
% (v,w could be vectors)
%
% the iqc's have the form
% v'*p*v>w'*p*w
% where p=x1+x2/(s+a(1))+x3/(s+a(2))+...>0,
% and also
% v'*q*w>0
% where
% q/2=z1+(y2*s+a(1)*z2)/(s^2-a(1)^2)+(y3*s+a(2)*z3)/(s^2-a(2)^2)+...
% xi arbitrary square matrices,
% yi arbitrary symmetric,
% zi arbitrary skew symmetric
%
% default a=1
%
% Written by ameg@mit.edu,  last modified October 13, 1997
% Last modified by cmj on 2013/5/4

global ABST

switch ABST.systemtype
    case 'continuous'
        if nargin<2
            a=1;
        end
        m=size(v,1);
        n=length(a);
        sys_s=tf([1 0],1);
        w=signal(m);
        if m>1,
            x{1}=rectangular(m);
            y{1}=OO(m);
            z{1}=skew(m);
            p=x{1};
            q=z{1};
            v'*(x{1}*v)>w'*(x{1}*w); %#ok<*VUNUS>
            v'*z{1}*w==0; %#ok<*EQEFF>
            for i=1:n,
                x{1+i}=rectangular(m); %#ok<*AGROW>
                y{1+i}=symmetric(m);
                z{1+i}=skew(m);
                Ga=abst(1/(sys_s+a(i)));
                va=Ga*v;
                wa=Ga*w;
                p=p+x{1+i}*Ga;
                q=q+y{1+i}*Ga-Ga'*y{1+i}+z{1+i}*Ga+Ga'*z{1+i};
                v'*(x{1+i}*va)>w'*x{1+i}*wa;
                v'*y{1+i}*wa>va'*(y{1+i}*w);
                v'*z{1+i}*wa+va'*z{1+i}*w>0;
            end
        else
            x{1}=rectangular;
            y{1}=OO;
            p=x{1};
            q=0;
            v'*(x{1}*v)>w'*(x{1}*w);
            for i=1:n,
                x{1+i}=rectangular;
                y{1+i}=symmetric;
                Ga=abst(1/(sys_s+a(i)));
                va=Ga*v;
                wa=Ga*w;
                p=p+x{1+i}*Ga;
                q=q+y{1+i}*Ga-Ga'*y{1+i};
                v'*(x{1+i}*va)>w'*x{1+i}*wa;
                v'*y{1+i}*wa==va'*(y{1+i}*w);
            end
        end
        p>0;
    case 'discrete'
        if nargin<2
            a=-0.3679;
        end
        m=size(v,1);
        n=length(a);
        sys_z=tf([1 0],1,-1);
        w=signal(m);
        if m>1,
            x{1}=rectangular(m);
            y{1}=OO(m);
            z{1}=skew(m);
            p=x{1};
            q=z{1};
            v'*(x{1}*v)>w'*(x{1}*w);
            v'*z{1}*w==0;
            for i=1:n,
                x{1+i}=rectangular(m);
                y{1+i}=symmetric(m);
                z{1+i}=skew(m);
                Ga=abst(1/(sys_z+a(i)));
                va=Ga*v;
                wa=Ga*w;
                p=p+x{1+i}*Ga;
                q=q+y{1+i}*Ga-Ga'*y{1+i}+z{1+i}*Ga+Ga'*z{1+i};
                v'*(x{1+i}*va)>w'*(x{1+i}*wa);
                v'*y{1+i}*wa>va'*(y{1+i}*w);
                v'*z{1+i}*wa+va'*z{1+i}*w>0;
            end
        else
            x{1}=rectangular;
            y{1}=OO;
            p=x{1};
            q=0;
            v'*(x{1}*v)>w'*x{1}*w;
            for i=1:n,
                x{1+i}=rectangular;
                y{1+i}=symmetric;
                Ga=abst(1/(sys_z+a(i)));
                va=Ga*v;
                wa=Ga*w;
                p=p+x{1+i}*Ga;
                q=q+y{1+i}*Ga-Ga'*y{1+i};
                v'*(x{1+i}*va)>w'*x{1+i}*wa;
                v'*y{1+i}*wa==va'*(y{1+i}*w);
            end
        end
        p>0;
end
