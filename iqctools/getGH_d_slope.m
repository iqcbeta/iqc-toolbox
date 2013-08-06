function [G,H]=getGH_d_slope(xxa,xxb,xxc,xxd,ddd,a,ign)
%
% USE
%   [G,H]=getGH_d_slope(xxa,xxb,xxc,xxd,ddd,a,ign);
%
% PURPOSE
%       Construct the multipliers of the iqc for repeated
%       slope-restricted nonlinearities.
% INPUTS
%    xxa: value of the variable "xa" (obtained by xxa=value_iqc(xa))
%    xxb: value of the variable "xb" (obtained by xxa=value_iqc(xb))
%    xxc: value of the variable "xc" (obtained by xxa=value_iqc(xc))
%    xxd: value of the variable "xd" (obtained by xxa=value_iqc(xd))
%    ddd: value of the variable "dd" (obtained by xxa=value_iqc(dd))
%         (All these variables are defined by the iqc_d_slope
%          function. It is assumed the command iqc_value was
%          issued after the optimization).
%      a: second argument used in iqc_d_slope
%    ign:fifth argument used in iqc_d_slope
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
% SEE ALSO: IQC_D_SLOPE, IQC_D_SLOPE_ODD, GETGH_D_SLOPE_ODD
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
countA=1;
countC=1;
countBD=1;

for n1=1:n
    for n2=n1:n
        m=length(a{n1,n2});
        for ndx=1:m
            aa=a{n1,n2}(ndx);
            keyboard
            
            if aa==Inf;
                G(n1,n2)=G(n1,n2)-xxa(countA);
                countA = countA+1;
            elseif isreal(aa)
                H(n1,n2)=H(n1,n2)+xxa(countA)*aa/(s+aa);
                countA=countA+1;
                H(n1,n2)=H(n1,n2)+xxc(countC)*aa/(-s+aa);
                countC=countC+1;
            elseif ~isreal(aa);
                b=abs(real(aa));
                c=abs(imag(aa));
                s1=b/(s+b);
                s2=b*c/((s+b)*(s+b)+c*c);
                s3=b*(s+b)/((s+b)*(s+b)+c*c);
                hh=s1+s2;
                m=s1-s2;
                r=s1+s3;
                t=s1-s3;
                s1_=b/(-s+b);
                s2_=b*c/((-s+b)*(-s+b)+c*c);
                s3_=b*(-s+b)/((-s+b)*(-s+b)+c*c);
                h_=s1_+s2_;
                m_=s1_-s2_;
                r_=s1_+s3_;
                t_=s1_-s3_;
                
                if ign==0
                    H(n1,n2) = H(n1,n2)+xxa(countA)*hh;
                    H(n1,n2) = H(n1,n2)+xxa(countA+1)*r;
                    countA = countA+2;
                    H(n1,n2) = H(n1,n2)+xxc(countC)*h_;
                    H(n1,n2) = H(n1,n2)+xxc(countC+1)*r_;
                    countC = countC+2;
                    H(n1,n2) = H(n1,n2)+xxb(countBD)*m;
                    H(n1,n2) = H(n1,n2)+xxb(countBD+1)*t;
                    H(n1,n2) = H(n1,n2)+xxd(countBD)*m_;
                    H(n1,n2) = H(n1,n2)+xxd(countBD+1)*t_;
                    countBD = countBD+2;
                    
                elseif ign==1;
                    H(n1,n2) = H(n1,n2)+xxa(countA)*r;
                    countA = countA+1;
                    H(n1,n2) = H(n1,n2)+xxc(countC)*r_;
                    countC = countC+1;
                    H(n1,n2) = H(n1,n2)+xxb(countBD)*t;
                    H(n1,n2) = H(n1,n2)+xxd(countBD)*t_;
                    countBD = countBD+1;
                    
                elseif ign==2;
                    H(n1,n2) = H(n1,n2)+xxa(countA)*hh;
                    countA = countA+1;
                    H(n1,n2) = H(n1,n2)+xxc(countC)*h_;
                    countC = countC+1;
                    H(n1,n2) = H(n1,n2)+xxb(countBD)*m;
                    H(n1,n2) = H(n1,n2)+xxd(countBD)*m_;
                    countBD = countBD+1;
                    
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
