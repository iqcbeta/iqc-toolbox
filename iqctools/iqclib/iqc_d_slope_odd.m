function [w,xa,xb,xc,xd,dd]=iqc_d_slope_odd(v,a,alpha,beta,ign)
% function [w,xa,xb,xc,xd,dd]=iqc_d_slope_odd(v,a,alpha,beta,ign)
%
% Defines iqc's for the relation w(t)=PHI(v(t)),
% where PHI is a DIAGONAL operator:
%
%  PHI= [phi 0......0;
%        0  phi 0...0;
%        :          :
%        0......0 phi];    and phi is odd and satisfies:
%
%           phi(x) - phi(y)
% alpha  < ---------------- < alpha+beta     for beta>0
%                x - y
%
% The IQC's have the form:
%
% \int q' D p  dt > \int q' H p dt
%
% where  p = (1+alpha/beta)v - 1/beta w   (auxiliary input)
%        q =  -alpha v +  w               (auxiliary output)
%
% D is diagonal and H=\H^T is a convolution operator
% with:     sum_j ||h_ij||_1<=D_ii,
%
% v : input signal
% a : SYMMETRIC (n x n) cell of vectors a_ij  containing
%     the (-)poles of h_ij
%
% IF a IS EITHER EMPTY OR NOT SPECIFIED, THEN H_ij=1 (static multipliers)
% IF a IS AN n x n  CELL, THEN
%  (1) if a{i,j}(k)==Inf            H_ij(k) = 1
%  (2) if a{i,j}(k)==a              H_ij(k) = a/(s + a)
%  (3) if a{i,j}(k)==b+i*c     i)  H_ij(k) = b*c/((s + b)^2+c^2)
%                              ii) H_ij(k+1) = b(s+b)/((s + b)^2+c^2)
%
% ign=1 : ignores multipliers (3)-i
% ign=2 : ignores multipliers (3)-ii
% ign=0 : uses both multipliers (3)-i and (3)-ii
%
% default:   H_ij=1, alpha=0,  beta=1 ign=0
%
% see also: IQC_D_SLOPE, IQC_SLOPE_ODD, IQC_SLOPE, IQC_MONOTONIC
%
% Written by fdamato@lids.mit.edu,    June 1998
% last modified Sept 19, 2000. July 24, 1998 by fdamato@ecn.purdue.edu
% Last modified by cmj on 2013/5/1

n=size(v,1);

if nargin<5
    ign=0;
end
if ign<0 || ign>2
    disp_str(69,'ign','0,1,2')
end
if nargin<4 || isempty(beta)
    beta=1;
end
if nargin<3 || isempty(alpha)
    alpha=0;
end
if nargin<2 || isempty(a)
    a=cell(n,n);
    for n1=1:n
        for n2=1:n
            a{n1,n2}=Inf;
        end
    end
else
    if size(a,1)~=n || size(a,2)~=n;
        error('inadequate dimensions of argument "a"');
    end
    for n1=1:n
        for n2=1:n
            if (isempty(a{n1,n2})+isempty(a{n2,n1}))==0;
                if sort(a{n1,n2})~=sort(a{n2,n1});
                    error('cell "a" must be symmetric');
                end
            elseif (isempty(a{n1,n2})+isempty(a{n2,n1}))==1;
                error('cell "a" must be symmetric');
            end
        end
    end
end
if nargin<1
    disp_str(1)
end

global ABST

switch ABST.systemtype
    case 'continuous'
        
        s = tf([1 0],1);
        w = signal(n);
        r=1+alpha/beta;
        p=r*v-1/beta*w;
        q=-alpha*v+w;
        
        Ninf=0;Npol=0;
        for n1=1:n;	% determine # of variables to define
            for n2=n1:n;
                aux=length(a{n1,n2});
                for ndx=1:aux;
                    if a{n1,n2}(ndx)==Inf;
                        Ninf=Ninf+1;
                    elseif isreal(a{n1,n2}(ndx)) || ign~=0;
                        Npol=Npol+1;
                    else 	% a complex and ign==0;
                        Npol=Npol+2;
                    end
                end % ndx
            end % n2
        end % n1
        
        dd=rectangular(1,n);		% variable definitions
        xa=rectangular(1,Ninf+Npol);
        xb=rectangular(1,Ninf+Npol);
        if Npol>0;
            xc=rectangular(1,Npol);
            xd=rectangular(1,Npol);
        end
        
        for ndx=1:Ninf+Npol;		% variables are positive
            xa{ndx}>0; %#ok<*VUNUS>
            xb{ndx}>0;
        end
        for ndx=1:Npol;
            xc{ndx}>0;
            xd{ndx}>0;
        end
        
        countAB=0;
        countCD=0;
        
        for n1=1:n
            n2=n1;
            i=num2str(n1);
            q(n1)'*dd(n1)*p(n1)>0;
            eval(['sumROW_' i '=dd(n1);']);
            
            m=length(a{n1,n1});          % -------diagonal  entries
            for ndx=1:m;
                aa=a{n1,n1}(ndx);
                countAB=countAB+1;
                if aa==Inf;		% static multiplier
                    q(n1)'*-xa{countAB}*p(n1)>0;
                    q(n1)'*xb{countAB}*p(n1)>0;
                    eval(['sumROW_' i '=sumROW_' i '-xa{countAB}-xb{countAB};']);
                elseif isreal(aa) % real pole
                    countCD=countCD+1;
                    Hij=aa/(s+aa);
                    HijPj=Hij*p(n1);
                    HijQj=Hij*q(n1);
                    
                    q(n1)'*-xa{countAB}*HijPj > 0;
                    q(n1)'*xb{countAB}*HijPj > 0;
                    p(n1)'*-xc{countCD}*HijQj > 0;
                    p(n1)'*xd{countCD}*HijQj > 0;
                    eval(['sumROW_' i '=sumROW_' i '-xa{countAB}-xb{countAB};']);
                    eval(['sumROW_' i '=sumROW_' i '-xc{countCD}-xd{countCD};']);
                else 	% complex poles
                    countCD=countCD+1;
                    b=real(aa);bb=abs(b);
                    c=imag(aa);cc=abs(c);
                    NormHij=bb*cc/(bb*bb+cc*cc)*(1+exp(-bb*pi/cc))/(1-exp(-bb*pi/cc));
                    if ign~=1;
                        Hij=b*c/((s+b)*(s+b)+c*c);
                        HijPj=Hij*p(n1);
                        HijQj=Hij*q(n1);
                        q(n1)'*-xa{countAB}*HijPj > 0;
                        q(n1)'*xb{countAB}*HijPj > 0;
                        p(n1)'*-xc{countCD}*HijQj > 0;
                        p(n1)'*xd{countCD}*HijQj > 0;
                        eval(['sumROW_' i '=sumROW_' i ...
                            '-(xa{countAB}+xb{countAB})*NormHij;']);
                        eval(['sumROW_' i '=sumROW_' i ...
                            '-(xc{countCD}+xd{countCD})*NormHij;']);
                    end
                    if ign==0;
                        countAB=countAB+1;
                        countCD=countCD+1;
                    end
                    if ign~=2;
                        Mij=b*(s+b)/((s+b)*(s+b)+c*c);
                        MijPj=Mij*p(n1);
                        MijQj=Mij*q(n1);
                        NormMij=NormHij*exp(-bb*pi/2/cc)+...
                            (bb*bb+cc*bb*exp(-bb*pi/2/cc))/(bb*bb+cc*cc);
                        q(n1)'*-xa{countAB}*MijPj > 0;
                        q(n1)'*xb{countAB}*MijPj > 0;
                        p(n1)'*-xc{countCD}*MijQj > 0;
                        p(n1)'*xd{countCD}*MijQj > 0;
                        eval(['sumROW_' i '=sumROW_' i ...
                            '-(xa{countAB}+xb{countAB})*NormMij;']);
                        eval(['sumROW_' i '=sumROW_' i ...
                            '-(xc{countCD}+xd{countCD})*NormMij;']);
                    end % ign==
                end % aa==Inf
            end % ndx
            
            for n2=n1+1:n       % ----------------off diagonal entries
                m=length(a{n1,n2});
                for ndx=1:m;
                    aa=a{n1,n2}(ndx);
                    countAB=countAB+1;
                    if aa==Inf;		% static multiplier
                        q(n1)'*-xa{countAB}*p(n2)>0;
                        q(n2)'*-xa{countAB}*p(n1)>0;
                        q(n1)'*xb{countAB}*p(n2)>0;
                        q(n2)'*xb{countAB}*p(n1)>0;
                        eval(['sumROW_' i '=sumROW_' i '-xa{countAB}-xb{countAB};']);
                    elseif isreal(aa) % real pole
                        countCD=countCD+1;
                        Hij=aa/(s+aa);
                        HijPj=Hij*p(n2);
                        HijQj=Hij*q(n2);
                        q(n1)'*-xa{countAB}*HijPj > 0;
                        q(n1)'*xb{countAB}*HijPj > 0;
                        p(n1)'*-xc{countCD}*HijQj > 0;
                        p(n1)'*xd{countCD}*HijQj > 0;
                        
                        HjiPi=Hij*p(n1);
                        HjiQi=Hij*q(n1);
                        q(n2)'*-xa{countAB}*HjiPi > 0;
                        q(n2)'*xb{countAB}*HjiPi > 0;
                        p(n2)'*-xc{countCD}*HjiQi > 0;
                        p(n2)'*xd{countCD}*HjiQi > 0;
                        eval(['sumROW_' i '=sumROW_' i '-xa{countAB}-xb{countAB};']);
                        eval(['sumROW_' i '=sumROW_' i '-xc{countCD}-xd{countCD};']);
                    else 	% complex poles
                        countCD=countCD+1;
                        b=real(aa);bb=abs(b);
                        c=imag(aa);cc=abs(c);
                        NormHij=bb*cc/(bb*bb+cc*cc)*(1+exp(-bb*pi/cc))/(1-exp(-bb*pi/cc));
                        if ign~=1;
                            Hij=b*c/((s+b)*(s+b)+c*c);
                            HijPj=Hij*p(n2);
                            HijQj=Hij*q(n2);
                            q(n1)'*-xa{countAB}*HijPj > 0;
                            q(n1)'*xb{countAB}*HijPj > 0;
                            p(n1)'*-xc{countCD}*HijQj > 0;
                            p(n1)'*xd{countCD}*HijQj > 0;
                            
                            HjiPi=Hij*p(n1);
                            HjiQi=Hij*q(n1);
                            q(n2)'*-xa{countAB}*HjiPi > 0;
                            q(n2)'*xb{countAB}*HjiPi > 0;
                            p(n2)'*-xc{countCD}*HjiQi > 0;
                            p(n2)'*xd{countCD}*HjiQi > 0;
                            
                            eval(['sumROW_' i '=sumROW_' i ...
                                '-(xa{countAB}+xb{countAB})*NormHij;']);
                            eval(['sumROW_' i '=sumROW_' i ...
                                '-(xc{countCD}+xd{countCD})*NormHij;']);
                        end
                        if ign==0;
                            countAB=countAB+1;
                            countCD=countCD+1;
                        end
                        if ign~=2;
                            Mij=b*(s+b)/((s+b)*(s+b)+c*c);
                            NormMij=NormHij*exp(-bb*pi/2/cc)+...
                                (bb*bb+cc*bb*exp(-bb*pi/2/cc))/(bb*bb+cc*cc);
                            MijPj=Mij*p(n2);
                            MijQj=Mij*q(n2);
                            q(n1)'*-xa{countAB}*MijPj > 0;
                            q(n1)'*xb{countAB}*MijPj > 0;
                            p(n1)'*-xc{countCD}*MijQj > 0;
                            p(n1)'*xd{countCD}*MijQj > 0;
                            
                            MjiPi=Mij*p(n1);
                            MjiQi=Mij*q(n1);
                            q(n2)'*-xa{countAB}*MjiPi > 0;
                            q(n2)'*xb{countAB}*MjiPi > 0;
                            p(n2)'*-xc{countCD}*MjiQi > 0;
                            p(n2)'*xd{countCD}*MjiQi > 0;
                            
                            eval(['sumROW_' i '=sumROW_' i ...
                                '-(xa{countAB}+xb{countAB})*NormMij;']);
                            eval(['sumROW_' i '=sumROW_' i ...
                                '-(xc{countCD}+xd{countCD})*NormMij;']);
                        end
                    end
                end % ndx
            end % n2
            
            for nfix1=1:n1-1  % contributions to sum_i from the lower triangular
                [countAB,countCD]=findindex(nfix1,n1,a,ign);
                aa = a{nfix1,n1};
                for ndx=1:length(aa);
                    countAB=countAB+1;
                    if aa(ndx)==Inf;
                        eval(['sumROW_' i '=sumROW_' i '-xa{countAB}-xb{countAB};']);
                    elseif isreal(aa(ndx))
                        countCD=countCD+1;
                        eval(['sumROW_' i '=sumROW_' i '-xa{countAB}-xb{countAB};']);
                        eval(['sumROW_' i '=sumROW_' i '-xc{countCD}-xd{countCD};']);
                    elseif ~isreal(aa(ndx))
                        countCD=countCD+1;
                        NormHij=bb*cc/(bb*bb+cc*cc)*(1+exp(-bb*pi/cc))/(1-exp(-bb*pi/cc));
                        NormMij=NormHij*exp(-bb*pi/2/cc)+...
                            (bb*bb+cc*bb*exp(-bb*pi/2/cc))/(bb*bb+cc*cc);
                        if ign==0
                            eval(['sumROW_' i '=sumROW_' i ...
                                '-(xa{countAB}+xb{countAB})*NormHij;']);
                            eval(['sumROW_' i '=sumROW_' i ...
                                '-(xc{countCD}+xd{countCD})*NormHij;']);
                            countAB=countAB+1;
                            countCD=countCD+1;
                            eval(['sumROW_' i '=sumROW_' i ...
                                '-(xa{countAB}+xb{countAB})*NormMij;']);
                            eval(['sumROW_' i '=sumROW_' i ...
                                '-(xc{countCD}+xd{countCD})*NormMij;']);
                        elseif ign==1
                            eval(['sumROW_' i '=sumROW_' i ...
                                '-(xa{countAB}+xb{countAB})*NormMij;']);
                            eval(['sumROW_' i '=sumROW_' i ...
                                '-(xc{countCD}+xd{countCD})*NormMij;']);
                        elseif ign==2
                            eval(['sumROW_' i '=sumROW_' i ...
                                '-(xa{countAB}+xb{countAB})*NormHij;']);
                            eval(['sumROW_' i '=sumROW_' i ...
                                '-(xc{countCD}+xd{countCD})*NormHij;']);
                        end % ign==0
                    end % if aa==Inf
                end % ndx
            end % nfix
        end % n1
        
        for n1=1:n;
            i=num2str(n1);
            eval(['sumROW_' i '>0;']);
        end
        
        if Npol==0; % assigns a(ny) value to xc, xd to avoid warnings
            xc=0;xd=0;
        end
        
    case 'discrete'
        disp_str(70,'iqc_d_slope_odd','discrete')
end

%% ----- auxiliay subroutine --------
function  [countAB,countCD]=findindex(ndx1,ndx2,a,ign)

if ndx1>ndx2;
    countAB=-1;
    countCD=-1;
    return
end
n=length(a);
countAB=0;
countCD=0;

for n2=1:n
    for n1=1:ndx1-1
        aa=a{n1,n2};
        m=length(aa);
        for ndx=1:m
            if aa(ndx)==Inf;
                countAB=countAB+1;
            elseif isreal(aa(ndx))
                countAB=countAB+1;
                countCD=countCD+1;
            elseif ~isreal(aa(ndx))
                if ign==0
                    countAB=countAB+2;
                    countCD=countCD+2;
                else
                    countAB=countAB+1;
                    countCD=countCD+1;
                end
            end
        end
    end
end

for n2=ndx1:ndx2-1
    aa=a{ndx1,n2};
    m=length(aa);
    for ndx=1:m
        if aa(ndx)==Inf;
            countAB=countAB+1;
        elseif isreal(aa(ndx))
            countAB=countAB+1;
            countCD=countCD+1;
        elseif ~isreal(aa(ndx))
            if ign==0
                countAB=countAB+2;
                countCD=countCD+2;
            else
                countAB=countAB+2;
                countCD=countCD+2;
            end
        end
    end
end
