function [w,Qr,Pr,Qi,Pi]=iqc_multi_harmonic(v,K,a)
% function [w,Qr,Pr,Qi,Pi]=iqc_multi_harmonic(v,K,a)
%
% Defines IQCs for operators defined as a multiplication with a harmonic
% multiplier in the time domain, i.e., the input output relation
%
% w(t)=Delta(t)*v(t)
%
% where
%                 [cos(w0*k1*t)I]
%                 [sin(w0*k1*t)I]
%        Delta(t)=[    :        ],
%                 [cos(w0*kN*t)I]
%                 [cos(w0*kN*t)I]
%
% and where 0<k1<k2<...<kN are integers, w0 is any real number, and I is the
% nxn identity matrix, where n is the dimension of v.
%
% The IQCs have the form
%
% [v]'   [v]
% [ ] *M*[ ]>0, where M=Real(U'*Q*U),
% [w]    [w]
%
% where
%
%        [I 0 0  0 0 ... 0  0]
%        [0 I iI 0 0 ... 0  0]
%      U=[        :          ]
%        [        :          ]
%        [0 0 0  0 0 ... I iI]
%
% and Q=Q' satisfies
%
% [A'*P*A-P  A'*P*B]
% [                ]+[C D]'*Q*[C D]>0
% [B'*P*A    B'*P*B]
%
% for some P=P'. Here (A,B,C,D) is a controllable realization of
%
%                 [   I   ]
%           1     [z^{k1}I]
% H(z)=----------=[   :   ]
%      (z+a)^{kN} [   :   ]
%                 [z^{kN}I]
%
% Furthermore, in order to use the parametrization tau*Delta_H we need
% to impose the constraints
%
%                          [M11  M12]
% M11>0 and M22<0, where M=[        ].
%                          [M12' M22]
%
% Inputs: K=[k1,k2,...,kN]
%         a>1
%
% Default: a=2, K=1
%
% See also iqc_multi_harmonic_gb
%
%
% Written by ulfj@mit.edu April 10, 1998
% Last modified by cmj on 2013/5/4

if nargin<2
    K=1;
end

global ABST

switch ABST.systemtype
    case 'continuous'
        if nargin<3
            a=2;
        end
        n=size(v,1);
        N=length(K);
        KN=K(N);
        I=eye(n);
        E=[I zeros(n,KN*n)];
        for j=1:N
            E=[E;zeros(n,n*K(j)) I zeros(n,n*(KN-K(j)))]; %#ok<*AGROW>
        end
        Ur=[I zeros(n,2*N*n)];
        Ui=zeros(n,(1+2*N)*n);
        for j=1:N
            Ur=[Ur;zeros(n,n*(1+2*(j-1))) I zeros(n,n*(2*(N-j)+1))];
            Ui=[Ui;zeros(n,n*(2+2*(j-1))) I zeros(n,n*2*(N-j))];
        end
        B=[zeros((KN-1)*n,n);I];
        D1=[zeros(KN*n,n);I];
        D=E*D1;
        as=[];
        for k=KN:-1:1
            as=[as -a^k*gamma(KN+1)/(gamma(k+1)*gamma(KN-k+1))*I];
        end
        A=[zeros((KN-1)*n,n) eye((KN-1)*n);as];
        C1=[eye(KN*n);as];
        C=E*C1;
        Qr=symmetric((N+1)*n);
        Qi=skew((N+1)*n);
        Pr=symmetric(KN*n);
        if KN*n==1
            Pi=0;
        else
            Pi=skew(KN*n);
        end
        [ A'*Pr*A-Pr A'*Pi*A-Pi  A'*Pr*B  A'*Pi*B;...
            -A'*Pi*A+Pi A'*Pr*A-Pr -A'*Pi*B  A'*Pr*B;...
            B'*Pr*A    B'*Pi*A     B'*Pr*B  B'*Pi*B;...
            -B'*Pi*A    B'*Pr*A    -B'*Pi*B  B'*Pr*B]+...
            [ C'*Qr*C    C'*Qi*C     C'*Qr*D  C'*Qi*D;...
            -C'*Qi*C    C'*Qr*C    -C'*Qi*D  C'*Qr*D;...
            D'*Qr*C    D'*Qi*C     D'*Qr*D  D'*Qi*D;...
            -D'*Qi*C    D'*Qr*C    -D'*Qi*D  D'*Qr*D]>0; %#ok<*VUNUS>
        M=Ur'*Qr*Ur+Ui'*Qr*Ui+Ui'*Qi*Ur-Ur'*Qi*Ui;
        %
        %To satisfy tau dependence
        %
        H1=[zeros(n,2*n*N);eye(2*n*N)];
        H1'*M*H1<0;
        H2=[ones(n);zeros(2*n*N,n)];
        H2'*M*H2>0;
        %
        %The IQC
        %
        w=signal(2*n*N);
        [v;w]'*M*[v;w]>0;
    case 'discrete'
        disp_str(70,'iqc_multi_harmonic','discrete')
end
