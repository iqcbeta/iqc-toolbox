function [gain,estimation]=...
    iqc_estimation_l2gain_CTrep_lmilab(w,v,y,p,q,weight)
% function [gain,estimation]=...
%     iqc_estimation_l2gain_CTrep_lmilab(w,v,y,p,q,weight)
%
%          ---------
%          |       |
%      ----| delta |----
%      |   |       |   |
%    p |   ---------   | q
%      |               |
%      |   ---------   |
%      ----|       |----
%          |       |         v     -----    -----  e
%  w ----> |   P   |---------------| - |----| W |---->
%          |       |-------        -----    -----
%          ---------   y  |          ^
%                         |  -----   |
%                         -->| E |----
%                            -----  hat{v}

global ABST
str={};

switch ABST.systemtype
    case 'discrete',
        disp_str(70,'iqc_estimation_l2gain_CTrep_lmilab','discrete')
end

lmiparameter=ABST.lmiparameter;

%vrb=0 to suppress the output
vrb=(lmiparameter(5)==0)||(lmiparameter(5)==777);

if vrb,
    if lmiparameter(5)>100,
        disp_str(50,'iqc_estimation_l2gain_CTrep_lmilab')
    else
        disp_str(51,'iqc_estimation_l2gain_CTrep_lmilab')
    end
end

%% first, process the iqc abst log
E=iqc_extract(w,v);
E=iqc_reduce(E);

%% Api [Bpiq Bpip Bpiw] Cpi [Dpiq Dpip Dpiw]
Psi=multiplier_ss(E,q,p,w);

Api=Psi.Api;
Bpiq=Psi.Bpi_q;
Bpip=Psi.Bpi_p;
Bpiw=Psi.Bpi_w;
Cpi=Psi.Cpi;
Dpiq=Psi.Dpi_q;
Dpip=Psi.Dpi_p;
Dpiw=Psi.Dpi_w;

Pi_pos=Psi.ss_pos;
p_pos=Psi.p_pos;
w_pos=Psi.w_pos;

%% A Bw Bp Cy Dyw Dyp Cv Dvw Dvp Cq Dqw Dqp

N_pos=setdiff(1:E.nstates,Pi_pos);

A=E.ab(N_pos,N_pos);
Bw=E.ab(N_pos,w_pos);
Bp=E.ab(N_pos,p_pos);

ny=double(y);
Cy=E.C{ny}(:,N_pos);
Dyw=E.C{ny}(:,w_pos);
Dyp=E.C{ny}(:,p_pos);

nv=double(v);
Cv=E.C{nv}(:,N_pos);
Dvw=E.C{nv}(:,w_pos);
Dvp=E.C{nv}(:,p_pos);

nq=double(q);
Cq=E.C{nq}(:,N_pos);
Dqw=E.C{nq}(:,w_pos);
Dqp=E.C{nq}(:,p_pos);

%%
Psi=ss(Api,[Bpiq Bpiw Bpip],Cpi,[Dpiq Dpiw Dpip]);

% Kalman decomposition
G1=ss(Psi.a,Psi.b(:,1:size(Cq,1)),Psi.c,Psi.d(:,1:size(Cq,1)));

[msys,u]=minreal(G1,[],false);

newA=u*(Psi.a)*u';
A1=newA(1:size(msys.a,1),1:size(msys.a,2));
A3=newA(1:size(msys.a,1),size(msys.a,2)+1:end);
A2=newA(size(msys.a,1)+1:end,size(msys.a,2)+1:end);

newB_1=u*Psi.b(:,1:size(Cq,1));
B1=newB_1(1:size(A1,1),:);

newB_2=u*Psi.b(:,size(Cq,1)+1:end);
B3=newB_2(1:size(A1,1),:);
B2=newB_2(size(A1,1)+1:end,:);

newC=(Psi.c)*u';
C1=newC(:,1:size(A1,2));
C2=newC(:,size(A1,2)+1:end);

D1=Psi.d(:,1:size(B1,2));
D2=Psi.d(:,size(B1,2)+1:end);

weight=ss(weight)*eye(size(Cv,1));
A_W=weight.a;
B_W=weight.b;
C_W=weight.c;
D_W=weight.d;

var=2;
lmi=4;

export_date=date;
str{1}=['%% Created by iqc_estimation_l2gain_CTrep_lmilab on ' export_date];
str{2}='\n%% Load coefficient matrices ...';
str{3}='load iqc_estimation_l2gain_CTrep_lmilab_exe';

%% initialize the LMI Control Toolbox
str{4}='\n%% Initialize the LMI Lab';
str{5}='setlmis([]);';
eval(str{5});

%% variable definite
defVar;

%% now we define the "non-KYP" lmi
nonKYP;

%% MpiR_*_* MpiL_*_*
multiplier_M;

%% Mpi11>0
k=E.nlmi+1;
ks=num2str(k);

Mpi11gt0;

%% now we define the "main" lmi
if vrb
    disp_str(54)
end
str{sc}='\n%% define the "main" lmi terms...'; %#ok<*NODEF>
sc=sc+1;

str{sc}='min_gamma=lmivar(1,[1 1]);';
eval(str{sc});
sc=sc+1;

na1=num2str(size(A_W,1));
str{sc}='hatA11=A_W;';
eval(str{sc});
sc=sc+1;

na2=num2str(size(A1,1));
na3=num2str(size(A2,1));
na4=num2str(size(A,1));
str{sc}=['hatA12=[zeros(',na1,',',cal_str(na2,na3),') B_W*Cv];'];
eval(str{sc});
sc=sc+1;
str{sc}=['hatA22=[A1 A3 B1*Cq;',...
    'zeros(',na3,',',na2,') A2 zeros(',na3,',',na4,');',...
    'zeros(',na4,',',cal_str(na2,na3),') A];'];
eval(str{sc});
sc=sc+1;

str{sc}=['hatA=[hatA11 hatA12;zeros(',cal_str(na2,na3,na4),...
    ',',na1,') hatA22];'];
eval(str{sc});
sc=sc+1;

nb1=num2str(size(Dvw,2));
nb2=num2str(size(Dvp,2));
str{sc}=['hatB=[B_W*Dvw B_W*Dvp;',...
    '[B1*Dqw B1*Dqp]+B3;',...
    'B2;',...
    'Bw Bp];'];
eval(str{sc});
sc=sc+1;
str{sc}=['hatB1=hatB(:,1:',nb1,');'];
eval(str{sc});
sc=sc+1;

str{sc}=['hatB2=hatB(:,',cal_str(nb1,'1'),':end);'];
eval(str{sc});
sc=sc+1;

nx1=num2str(size(B_W,2));
ny2=num2str(size(hatA22,2));
str{sc}=['hatB=[B_W;zeros(',ny2,',',nx1,')];'];
eval(str{sc});
sc=sc+1;

nc12=num2str(size(C1,1));
str{sc}=['hatC12=[zeros(',nc12,',',na1,') C1 C2 D1*Cq];'];
eval(str{sc});
sc=sc+1;

nc3=num2str(size(C_W,1));
str{sc}=['hatC3=[C_W zeros(',nc3,',',cal_str(na2,na3),') D_W*Cv];'];
eval(str{sc});
sc=sc+1;

str{sc}='hatE2=D_W;';
eval(str{sc});
sc=sc+1;

str{sc}='hatD12=[D1*Dqw D1*Dqp]+D2;';
eval(str{sc});
sc=sc+1;

str{sc}='hatD3=[D_W*Dvw D_W*Dvp];';
eval(str{sc});
sc=sc+1;

nx1=num2str(size(Cy,1));
str{sc}=['hatC=[zeros(',nx1,',',cal_str(na1,na2,na3),') Cy];'];
eval(str{sc});
sc=sc+1;

str{sc}='hatF1=Dyw;';
eval(str{sc});
sc=sc+1;

str{sc}='hatF2=Dyp;';
eval(str{sc});
sc=sc+1;

str{sc}=['[X,nX,sX]=lmivar(1,[',num2str(size(hatA,1)),' 1]);'];
eval(str{sc});
sc=sc+1;

nx1=num2str(size(hatA11,1));
nx2=num2str(size(hatA22,1));

if size(hatA11,1)~=0
    
    str{sc}=['[T11,nT11,sT11]=lmivar(1,[',nx1,' 1]);'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['[T12,nT12,sT12]=lmivar(2,[',nx1,' ',nx2,']);'];
    eval(str{sc});
    sc=sc+1;
    
end

str{sc}=['[T22,nT22,sT22]=lmivar(1,[',nx2,' 1]);'];
eval(str{sc});
sc=sc+1;

%%

% define hatM hatK L N
nx1=num2str(size(hatA,1));
str{sc}=['[hatK,nhatK,shatK]=lmivar(2,[',nx1,' ',nx1,']);'];
eval(str{sc});
sc=sc+1;

nx2=num2str(size(Cy,1));
str{sc}=['[L,nL,sL]=lmivar(2,[',nx1,' ',nx2,']);'];
eval(str{sc});
sc=sc+1;

ny2=num2str(size(Cv,1));
str{sc}=['[hatM,nhatM,shatM]=lmivar(2,[',ny2,' ',nx1,']);'];
eval(str{sc});
sc=sc+1;

str{sc}=['[N,nN,sN]=lmivar(2,[',ny2,',',nx2,']);'];
eval(str{sc});
sc=sc+1;

nx1=num2str(size(hatA11,1));
nx2=num2str(size(hatA22,1));

ouf_str1=['[eye(',nx1,') zeros(',nx1,',',nx2,')]'];
ouf_str2=['[zeros(',nx2,',',nx1,') eye(',nx2,')]'];

% Ycl>0
% [T11   0  | I    T12 ]
% [ 0   T22 | 0    T22']
% [--------------------]  > 0
% [ I    0  | X11  X12 ]
% [T12' T22 | X12' X22 ]

ks=num2str(str2double(ks)+1);

if size(hatA11,1)~=0
    
    str{sc}=['lmiterm([-',ks,' 1 1 T11],',ouf_str1,''',',ouf_str1,');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([-',ks,' 1 2 T12],',ouf_str1,''',',ouf_str2,');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([-',ks,' 1 2 0],',ouf_str1,'''*',ouf_str1,');'];
    eval(str{sc});
    sc=sc+1;
    
end

str{sc}=['lmiterm([-',ks,' 1 1 T22],',ouf_str2,''',',ouf_str2,');'];
eval(str{sc});
sc=sc+1;

str{sc}=['lmiterm([-',ks,' 1 2 -T22],',ouf_str2,''',',ouf_str2,');'];
eval(str{sc});
sc=sc+1;

str{sc}=['lmiterm([-',ks,' 2 2 X],1,1);'];
eval(str{sc});
sc=sc+1;

%% main lmi
% [[(Ae+Apsi)+(Ae+Apsi)'+C'*Pr0*C   Be+Bpsi+C'*Pr0*D]   *  ]
% [[                *                   D'*Pr0*D    ]   *  ] < 0
% [                   Theta11*[C3 D3]                   -I ]

ks=num2str(str2double(ks)+1);

str{sc}=['Theta11=eye(',nc3,');'];
eval(str{sc});
sc=sc+1;

str{sc}=['Theta12=zeros(',nc3,',',nb1,');'];
eval(str{sc});
sc=sc+1;

%% Ae+Ae'
% T1'*hatA*T2
if size(hatA11,1)~=0
    str{sc}=['lmiterm([',ks,' 1 1 T11],',ouf_str1,...
        '''*hatA11,',ouf_str1,',''s'');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 1 1 0],',ouf_str1,...
        '''*hatA12*',ouf_str2,'+(',ouf_str1,...
        '''*hatA12*',ouf_str2,')'');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 1 1 T12],-',ouf_str1,...
        '''*hatA11,',ouf_str2,',''s'');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 1 1 T12],',ouf_str1,...
        ''',hatA22*',ouf_str2,',''s'');'];
    eval(str{sc});
    sc=sc+1;
    
end

str{sc}=['lmiterm([',ks,' 1 1 T22],',ouf_str2,...
    ''',hatA22*',ouf_str2,',''s'');'];
eval(str{sc});
sc=sc+1;

% hatB*hatM
str{sc}=['lmiterm([',ks,' 1 1 hatM],hatB,1,''s'');'];
eval(str{sc});
sc=sc+1;

% T1'*hatA
if size(hatA11,1)~=0
    str{sc}=['lmiterm([',ks,' 1 2 0],',ouf_str1,...
        '''*',ouf_str1,'*hatA);'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 1 2 T12],',ouf_str1,...
        ''',',ouf_str2,'*hatA);'];
    eval(str{sc});
    sc=sc+1;
end
str{sc}=['lmiterm([',ks,' 1 2 T22],',ouf_str2,...
    ''',',ouf_str2,'*hatA);'];
eval(str{sc});
sc=sc+1;

% hatB*N*hatC
str{sc}=['lmiterm([',ks,' 1 2 N],hatB,hatC);'];
eval(str{sc});
sc=sc+1;

% hatK'
str{sc}=['lmiterm([',ks,' 1 2 -hatK],1,1);'];
eval(str{sc});
sc=sc+1;

% X*hatA
str{sc}=['lmiterm([',ks,' 2 2 X],1,hatA,''s'');'];
eval(str{sc});
sc=sc+1;

% L*hatC
str{sc}=['lmiterm([',ks,' 2 2 L],1,hatC,''s'');'];
eval(str{sc});
sc=sc+1;

%% Apsi+Apsi'
ouf_str3=['[zeros(',na2,',',na1,') eye(',na2,') zeros(',na2,',',...
    cal_str(na3,na4),')]'];

if size(A1,1)~=0
    str{sc}=['lmiterm([',ks,' 1 1 Xpsi],',ouf_str3,...
        ''',',ouf_str3,'*hatA,''s'');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 1 2 Xpsi],',ouf_str3,...
        ''',',ouf_str3,'*hatA,''s'');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 2 2 Xpsi],',ouf_str3,...
        ''',',ouf_str3,'*hatA,''s'');'];
    eval(str{sc});
    sc=sc+1;
end

%% Be+Bpsi
% T1'*[hatB1 hatB2]
if size(hatA11,1)~=0
    str{sc}=['lmiterm([',ks,' 1 3 0],',ouf_str1,...
        '''*',ouf_str1,'*[hatB1 hatB2]);'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 1 3 T12],',ouf_str1,...
        ''',',ouf_str2,'*[hatB1 hatB2]);'];
    eval(str{sc});
    sc=sc+1;
end
str{sc}=['lmiterm([',ks,' 1 3 T22],',ouf_str2,...
    ''',',ouf_str2,'*[hatB1 hatB2]);'];
eval(str{sc});
sc=sc+1;

% hatB*N*[hatF1 hatF2]
str{sc}=['lmiterm([',ks,' 1 3 N],hatB,[hatF1 hatF2]);'];
eval(str{sc});
sc=sc+1;

% X*[hatB1 hatB2]
str{sc}=['lmiterm([',ks,' 2 3 X],1,[hatB1 hatB2]);'];
eval(str{sc});
sc=sc+1;

% L*[hatF1 hatF2]
str{sc}=['lmiterm([',ks,' 2 3 L],1,[hatF1 hatF2]);'];
eval(str{sc});
sc=sc+1;

% Bpsi
if size(A1,1)~=0
    str{sc}=['lmiterm([',ks,' 1 3 Xpsi],',ouf_str3,...
        ''',',ouf_str3,'*[hatB1 hatB2]);'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 2 3 Xpsi],',ouf_str3,...
        ''',',ouf_str3,'*[hatB1 hatB2]);'];
    eval(str{sc});
    sc=sc+1;
end

%% C'*Pr0*D

ouf_str4=['[eye(',nb1,') zeros(',nb1,',',nb2,')]'];

% T2'*hatC3'*Theta12*[I 0]
if size(hatA11,1)~=0
    str{sc}=['lmiterm([',ks,' 1 3 T11],',ouf_str1,...
        ''',',ouf_str1,'*hatC3''*Theta12*',ouf_str4,');'];
    eval(str{sc});
    sc=sc+1;
    str{sc}=['lmiterm([',ks,' 1 3 -T12],-',ouf_str2,...
        ''',',ouf_str1,'*hatC3''*Theta12*',ouf_str4,');'];
    eval(str{sc});
    sc=sc+1;
end
str{sc}=['lmiterm([',ks,' 1 3 0],',ouf_str2,...
    '''*',ouf_str2,'*hatC3''*Theta12*',ouf_str4,');'];
eval(str{sc});
sc=sc+1;

% hatM'*hatE2'
str{sc}=['lmiterm([',ks,' 1 3 -hatM],1,hatE2''*Theta12*',ouf_str4,');'];
eval(str{sc});
sc=sc+1;

% hatC3'*Theta12*[I 0]
str{sc}=['lmiterm([',ks,' 2 3 0],hatC3''*Theta12*',ouf_str4,');'];
eval(str{sc});
sc=sc+1;

% hatC'*N'*hatE2'*Theta12*[I 0]
str{sc}=['lmiterm([',ks,' 2 3 -N],hatC'',hatE2''*Theta12*',ouf_str4,');'];
eval(str{sc});
sc=sc+1;

%% D'*Pr0*D
% [I;0]*Theta12'*D3
str{sc}=['lmiterm([',ks,' 3 3 0],(',ouf_str4,...
    '''*Theta12''*hatD3)+(',ouf_str4,...
    '''*Theta12''*hatD3)'');'];
eval(str{sc});
sc=sc+1;

str{sc}=['lmiterm([',ks,' 3 3 N],',ouf_str4,...
    '''*Theta12''*hatE2,[hatF1 hatF2],''s'');'];
eval(str{sc});
sc=sc+1;

% [I;0]*Theta22*[I 0]

str{sc}=['lmiterm([',ks,' 3 3 min_gamma],-',ouf_str4,''',',ouf_str4,');'];
eval(str{sc});
sc=sc+1;

%% Theta11^(1/2)*C3

if size(hatA11,1)~=0
    str{sc}=['lmiterm([',ks,' 4 1 T11],sqrtm(Theta11)*hatC3*',ouf_str1,...
        ''',',ouf_str1,');'];
    eval(str{sc});
    sc=sc+1;
    str{sc}=['lmiterm([',ks,' 4 1 T12],-sqrtm(Theta11)*hatC3*',ouf_str1,...
        ''',',ouf_str2,');'];
    eval(str{sc});
    sc=sc+1;
end
str{sc}=['lmiterm([',ks,' 4 1 0],sqrtm(Theta11)*hatC3*',ouf_str2,...
    '''*',ouf_str2,');'];
eval(str{sc});
sc=sc+1;

str{sc}=['lmiterm([',ks,' 4 1 hatM],sqrtm(Theta11)*hatE2,1);'];
eval(str{sc});
sc=sc+1;

str{sc}=['lmiterm([',ks,' 4 2 0],sqrtm(Theta11)*hatC3);'];
eval(str{sc});
sc=sc+1;

str{sc}=['lmiterm([',ks,' 4 2 N],sqrtm(Theta11)*hatE2,hatC);'];
eval(str{sc});
sc=sc+1;

%% Theta11^(1/2)*D3
str{sc}=['lmiterm([',ks,' 4 3 0],sqrtm(Theta11)*hatD3);'];
eval(str{sc});
sc=sc+1;

str{sc}=['lmiterm([',ks,' 4 3 N],sqrtm(Theta11)*hatE2,[hatF1 hatF2]);'];
eval(str{sc});
sc=sc+1;

%% -I
str{sc}=['lmiterm([',ks,' 4 4 0],-eye(',nc3,'));'];
eval(str{sc});
sc=sc+1;

%% (*)'*Pr0*[C D]
for j=1:E.nsimple,
    js=num2str(j);
    nh=0;
    ng=0;
    n=E.simples(j);
    for i=1:size(E.X{n},2),
        is=num2str(i);
        vrs = num2str(E.X{n}(1,i));
        
        % (*)'Mpi*[C1;C2]
        str{sc}=['lmiterm([',ks,' 1 1 ',vrs,...
            '],hatC12''*MpiL_',js,'_',is,...
            ',MpiR_',js,'_',is,'*hatC12,''s'');']; %#ok<*AGROW>
        eval(str{sc})
        sc=sc+1;
        
        str{sc}=['lmiterm([',ks,' 1 2 ',vrs,...
            '],hatC12''*MpiL_',js,'_',is,...
            ',MpiR_',js,'_',is,'*hatC12);'];
        eval(str{sc})
        sc=sc+1;
        
        str{sc}=['lmiterm([',ks,' 1 2 -',vrs,...
            '],hatC12''*MpiR_',js,'_',is,...
            ''',MpiL_',js,'_',is,'''*hatC12);'];
        eval(str{sc})
        sc=sc+1;
        
        str{sc}=['lmiterm([',ks,' 2 2 ',vrs,...
            '],hatC12''*MpiL_',js,'_',is,...
            ',MpiR_',js,'_',is,'*hatC12,''s'');'];
        eval(str{sc})
        sc=sc+1;
        
        % [C1;C2]'*Mpi*[D1;D2]
        str{sc}=['lmiterm([',ks,' 1 3 ',vrs,...
            '],hatC12''*MpiL_',js,'_',is,...
            ',MpiR_',js,'_',is,'*hatD12);'];
        eval(str{sc})
        sc=sc+1;
        
        str{sc}=['lmiterm([',ks,' 1 3 -',vrs,...
            '],hatC12''*MpiR_',js,'_',is,...
            ''',MpiL_',js,'_',is,'''*hatD12);'];
        eval(str{sc})
        sc=sc+1;
        
        str{sc}=['lmiterm([',ks,' 2 3 ',vrs,...
            '],hatC12''*MpiL_',js,'_',is,...
            ',MpiR_',js,'_',is,'*hatD12);'];
        eval(str{sc})
        sc=sc+1;
        
        str{sc}=['lmiterm([',ks,' 2 3 -',vrs,...
            '],hatC12''*MpiR_',js,'_',is,...
            ''',MpiL_',js,'_',is,'''*hatD12);'];
        eval(str{sc})
        sc=sc+1;
        
        % [D1;D2]'*Mpi*[D1;D2]
        str{sc}=['lmiterm([',ks,' 3 3 ',vrs,...
            '],hatD12''*MpiL_',js,'_',is,...
            ',MpiR_',js,'_',is,'*hatD12,''s'');'];
        eval(str{sc})
        sc=sc+1;
        
        %         testf
        
        nh=nh+E.X{n}(3,i);
        ng=ng+E.X{n}(2,i);
    end
end

%% solving the system of LMI's
str{sc}='\n%% Solving the system of LMIs ...';
sc=sc+1;
str{sc}='lmi=getlmis;';
eval(str{sc});
sc=sc+1;
str{sc}='ndec=decnbr(lmi);';
eval(str{sc});
sc=sc+1;

if vrb
    disp_str(56,num2str(ndec))
end

str{sc}='c=zeros(ndec,1);';
eval(str{sc});
sc=sc+1;
str{sc}='nc=decinfo(lmi,min_gamma);';
eval(str{sc});
sc=sc+1;
str{sc}='c(nc)=1;';
eval(str{sc});
sc=sc+1;
gain = inf;
str{sc}='[gain,xopt]=mincx(lmi,c,lmiparameter);';
eval(str{sc});
sc=sc+1;
if ~isempty(gain) && gain > 0
    if vrb
        disp_str(57,'L2-gain','Feasible')
    end
    str{sc}='gain=sqrt(gain)';
    eval([str{sc},';']);
    sc=sc+1;
else
    if vrb
        disp_str(57,'L2-gain','Infeasible')
    end
    str{sc}='gain=inf';
    eval([str{sc},';']);
    sc=sc+1;
end

estimation=[];
if ~isinf(gain)
    str{sc}='hatK=decs2mats(xopt,shatK);';
    eval(str{sc})
    sc=sc+1;
    str{sc}='hatM=decs2mats(xopt,shatM);';
    eval(str{sc})
    sc=sc+1;
    str{sc}='N=decs2mats(xopt,sN);';
    eval(str{sc})
    sc=sc+1;
    str{sc}='L=decs2mats(xopt,sL);';
    eval(str{sc})
    sc=sc+1;
    if size(hatA11,1)~=0
        str{sc}='T11=decs2mats(xopt,sT11);';
        eval(str{sc})
        sc=sc+1;
        str{sc}='T12=decs2mats(xopt,sT12);';
        eval(str{sc})
        sc=sc+1;
    else
        str{sc}='T11=zeros(0,0);';
        eval(str{sc})
        sc=sc+1;
        str{sc}=['T12=zeros(0,',num2str(size(hatA12,2)),');'];
        eval(str{sc})
        sc=sc+1;
    end
    str{sc}='T22=decs2mats(xopt,sT22);';
    eval(str{sc})
    sc=sc+1;
    str{sc}='X=decs2mats(xopt,sX);';
    eval(str{sc})
    sc=sc+1;
    
    nx1=num2str(size(T11,1));
    nx2=num2str(size(T22,1));
    str{sc}=['T1=[eye(',nx1,') zeros(',nx1,',',nx2,');T12'' T22];'];
    eval(str{sc})
    sc=sc+1;
    
    getEST;
    if any(eig(AE)>0)
        disp_str(59,'estimation unstable')
        return
    end
    
    str{sc}='estimation=ss(AE,BE,CE,DE);';
    eval(str{sc});
    sc=sc+1;
end

if lmiparameter(5)>100;
    save iqc_estimation_l2gain_CTrep_lmilab_exe A A1 A2 A3 B1 B2 B3...
        Dqp Dqw Bp Bw C1 C2 D1 D2 Cq Cv Cy Dvp Dvw Dyp Dyw...
        lmiparameter A_W B_W C_W D_W
    if exist('state1','var')
        save iqc_estimation_l2gain_CTrep_lmilab_exe state* -append
    end
    if exist('struct1','var')
        save iqc_estimation_l2gain_CTrep_lmilab_exe struct* -append
    end
    if exist('l_1_1','var')
        save iqc_estimation_l2gain_CTrep_lmilab_exe l_* -append
    end
    if exist('r_1_1','var')
        save iqc_estimation_l2gain_CTrep_lmilab_exe r_* -append
    end
    if exist('MpiL_1_1','var')
        save iqc_estimation_l2gain_CTrep_lmilab_exe MpiL_* -append
    end
    if exist('MpiR_1_1','var')
        save iqc_estimation_l2gain_CTrep_lmilab_exe MpiR_* -append
    end
    fid=fopen('iqc_estimation_l2gain_CTrep_lmilab_exe.m','wt');
    for i=1:sc-1
        fprintf(fid,[str{i} '\n']);
    end
    fclose(fid);
end