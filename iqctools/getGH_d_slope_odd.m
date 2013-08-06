function [G,H]=getGH_d_slope_odd(xxa,xxb,xxc,xxd,ddd,a,ign)
%
% USE
%   [G,H]=getGH_d_slope_odd(xxa,xxb,xxc,xxd,ddd,a,ign);
%
% PURPOSE
%       Construct the multipliers of the iqc for repeated
%       slope-restricted odd nonlinearities.
%
% INPUTS
%    xxa: value of the variable "xa" (obtained by xxa=value_iqc(xa))
%    xxb: value of the variable "xb" (obtained by xxa=value_iqc(xa))
%    xxc: value of the variable "xc" (obtained by xxa=value_iqc(xa))
%    xxd: value of the variable "xd" (obtained by xxa=value_iqc(xa))
%    ddd: value of the variable "dd" (obtained by xxa=value_iqc(xa))
%         (All these variables are defined by the iqc_d_slope_odd
%          function. It is assumed the command iqc_value was issued
%          after the optimization).
%      a: second argument used in iqc_d_slope_odd
%    ign: fifth argument used in iqc_d_slope_odd
%
% OUTPUTS
%   G: static multiplier, square n x n symmetric matrix
%   H: dynamic multiplier, square n x n symmetric transfer function
%
%   These multipliers are defined in Theorem 1 of
%   D'Amato, Rotea, Megretski, Jonsson, "new resluts for
%   analysis of systems with repeated nonlinearities"
%   Automatica, to appear.
%
% SEE ALSO: IQC_D_SLOPE, IQC_D_SLOPE_ODD, GETGH_D_SLOPE
%
%
% Written by fdamato@ecn.purdue.edu,    Sept 2000

%
% change call of "value" to "value_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu



if nargin==6;
    ign=0;
end

n=length(a);
G=diag(ddd);
H=tf(ss([],[],[],zeros(n,n)));

s=tf([1 0],1);
countAB=1;
countCD=1;

for n1=1:n
    for n2=n1:n
        m=length(a{n1,n2});
        for ndx=1:m
            aa=a{n1,n2}(ndx);
            if aa==Inf;
                G(n1,n2)=G(n1,n2)-xxa(countAB)+xxb(countAB);
                countAB = countAB+1;
            elseif isreal(aa)
                H(n1,n2)=H(n1,n2)+(xxa(countAB)-xxb(countAB))*aa/(s+aa);
                countAB=countAB+1;
                H(n1,n2)=H(n1,n2)+(xxc(countCD)-xxd(countCD))*aa/(-s+aa);
                countCD=countCD+1;
            elseif ~isreal(aa);
                b=real(aa);
                c=imag(aa);
                if ign==0
                    h1 = (b*c)/((s+b)^2+c^2);
                    h2 = (b*(s+b))/((s+b)^2+c^2);
                    H(n1,n2) = H(n1,n2)+(xxa(countAB)-xxb(countAB))*h1;
                    H(n1,n2) = H(n1,n2)+(xxa(countAB+1)-xxb(countAB+1))*h2;
                    countAB = countAB+2;
                    h3 = (b*c)/((-s+b)^2+c^2);
                    h4 = (b*(-s+b))/((-s+b)^2+c^2);
                    H(n1,n2) = H(n1,n2)+(xxc(countCD)-xxd(countCD))*h3;
                    H(n1,n2) = H(n1,n2)+(xxc(countCD+1)-xxd(countCD+1))*h4;
                    countCD = countCD+2;
                elseif ign==1;
                    h2 = (b*(s+b))/((s+b)^2+c^2);
                    H(n1,n2) = H(n1,n2)+(xxa(countAB)-xxb(countAB))*h2;
                    countAB = countAB+1;
                    h4 = (b*(-s+b))/((-s+b)^2+c^2);
                    H(n1,n2) = H(n1,n2)+(xxc(countCD)-xxd(countCD))*h4;
                    countCD = countCD+1;
                elseif ign==2;
                    h1 = (b*c)/((s+b)^2+c^2);
                    H(n1,n2) = H(n1,n2)+(xxa(countAB)-xxb(countAB))*h1;
                    countAB = countAB+1;
                    h3 = (b*c)/((-s+b)^2+c^2);
                    H(n1,n2) = H(n1,n2)+(xxc(countCD)-xxd(countCD))*h3;
                    countCD = countCD+1;
                end
            end
        end
    end
end

G=G+transpose(G)-diag(diag(G));

for ndx=2:n
    for n2=1:ndx-1;
        H(ndx,n2)=H(n2,ndx);
    end
end
