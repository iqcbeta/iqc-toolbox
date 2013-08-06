function [gain,estimation]=...
    iqc_estimation_l2gain_DTeli_lmilab(w,v,y,p,q,weight)
% function [gain,estimation]=...
%     iqc_estimation_l2gain_DTeli_lmilab(w,v,y,p,q,weight)
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
    case 'continuous',
        disp_str(70,'iqc_estimation_l2gain_DTeli_lmilab','continuous')
end

lmiparameter=ABST.lmiparameter;

%vrb=0 to suppress the output
vrb=(lmiparameter(5)==0)||(lmiparameter(5)==777);

if vrb,
    if lmiparameter(5)>100,
        disp_str(50,'iqc_estimation_l2gain_DTeli_lmilab')
    else
        disp_str(51,'iqc_estimation_l2gain_DTeli_lmilab')
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
Psi=ss(Api,[Bpiq Bpiw Bpip],Cpi,[Dpiq Dpiw Dpip],-1);

% Kalman decomposition
G1=ss(Psi.a,Psi.b(:,1:size(Cq,1)),Psi.c,Psi.d(:,1:size(Cq,1)),-1);

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
str{1}=['%% Created by iqc_estimation_l2gain_DTeli_lmilab on ' export_date];
str{2}='\n%% Load coefficient matrices ...';
str{3}='load iqc_estimation_l2gain_DTeli_lmilab_exe';

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

ks=num2str(str2double(ks)+1);

str{sc}=['Theta11=eye(',nc3,');'];
eval(str{sc});
sc=sc+1;

str{sc}=['Theta12=zeros(',nc3,',',nb1,');'];
eval(str{sc});
sc=sc+1;

%% right lmi
str{sc}='\n%% define the "main" right lmi terms...';
sc=sc+1;

nx1=num2str(size(hatA,2));
nx2=num2str(size(hatB1,2));
nx3=num2str(size(hatB2,2));
str{sc}=['ouf_main_right=[eye(',nx1,') zeros(',nx1,',',...
    cal_str(nx2,nx3),');',...
    'hatA hatB1 hatB2;',...
    'hatC12 hatD12;',...
    'hatC3 hatD3;',...
    'zeros(',nx2,',',nx1,') eye(',nx2,') zeros(',nx2,',',nx3,')];'];
eval(str{sc});
sc=sc+1;

str{sc}='NWr=[hatC hatF1 hatF2];';
eval(str{sc});
sc=sc+1;

str{sc}='Wr=null(NWr,''r'');';
eval(str{sc});
sc=sc+1;

ouf_str3=['[zeros(',na2,',',na1,') eye(',na2,') zeros(',na2,',',...
    cal_str(na3,na4),')]'];

str{sc}=['lmiterm([',ks,' 1 1 X],-1,1);'];
eval(str{sc});
sc=sc+1;

if size(A1,1)~=0
    str{sc}=['lmiterm([',ks,' 1 1 Xpsi],-',ouf_str3,''',',ouf_str3,');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 2 2 Xpsi],',ouf_str3,''',',ouf_str3,');'];
    eval(str{sc});
    sc=sc+1;
end

str{sc}=['lmiterm([',ks,' 2 2 X],1,1);'];
eval(str{sc});
sc=sc+1;



str{sc}=['lmiterm([',ks,' 4 4 0],Theta11);'];
eval(str{sc});
sc=sc+1;

str{sc}=['lmiterm([',ks,' 4 5 0],Theta12);'];
eval(str{sc});
sc=sc+1;

str{sc}=['lmiterm([',ks,' 5 5 min_gamma],-1,eye(',nb1,'));'];
eval(str{sc});
sc=sc+1;

for j=1:E.nsimple,
    js=num2str(j);
    nh=0;
    ng=0;
    n=E.simples(j);
    for i=1:size(E.X{n},2),
        is=num2str(i);
        vrs = num2str(E.X{n}(1,i));
        
        % (*)'Mpi*[C1;C2]
        str{sc}=['lmiterm([',ks,' 3 3 ',vrs,...
            '],MpiL_',js,'_',is,...
            ',MpiR_',js,'_',is,',''s'');'];
        eval(str{sc})
        sc=sc+1;
        
        nh=nh+E.X{n}(3,i);
        ng=ng+E.X{n}(2,i);
    end
end

str{sc}=['lmiterm([',ks,' 0 0 0],ouf_main_right*Wr);'];
eval(str{sc});
sc=sc+1;

%% left lmi

ks=num2str(str2double(ks)+1);

%% -(bY+bXpsi)+A1'*Xpsi*A1
% -bY
if size(hatA11,1)~=0
    str{sc}=['lmiterm([',ks,' 1 1 T11],-',ouf_str1,...
        ''',',ouf_str1,');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 2 1 0],-',ouf_str1,...
        ''',',ouf_str1,');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 2 1 -T12],-',ouf_str2,...
        ''',',ouf_str1,');'];
    eval(str{sc});
    sc=sc+1;
end
str{sc}=['lmiterm([',ks,' 1 1 T22],-',ouf_str2,...
    ''',',ouf_str2,');'];
eval(str{sc});
sc=sc+1;

str{sc}=['lmiterm([',ks,' 2 1 T22],-',ouf_str2,...
    ''',',ouf_str2,');'];
eval(str{sc});
sc=sc+1;

str{sc}=['lmiterm([',ks,' 2 2 X],-1,1);'];
eval(str{sc});
sc=sc+1;

% bXpsi
ouf_str3=['[zeros(',na2,',',na1,') eye(',na2,') zeros(',na2,',',...
    cal_str(na3,na4),')]'];

if size(A1,1)~=0
    str{sc}=['lmiterm([',ks,' 1 1 Xpsi],-',ouf_str3,...
        ''',',ouf_str3,');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 2 1 Xpsi],-',ouf_str3,...
        ''',',ouf_str3,');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 2 2 Xpsi],-',ouf_str3,...
        ''',',ouf_str3,');'];
    eval(str{sc});
    sc=sc+1;
end

% A1'*Xpsi*A1
ouf_str5=['[zeros(',na2,',',na1,') A1 A3 B1*Cq]'];

if size(A1,1)~=0
    str{sc}=['lmiterm([',ks,' 1 1 Xpsi],',ouf_str5,...
        ''',',ouf_str5,');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 2 1 Xpsi],',ouf_str5,...
        ''',',ouf_str5,');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 2 2 Xpsi],',ouf_str5,...
        ''',',ouf_str5,');'];
    eval(str{sc});
    sc=sc+1;
end

%% A1'*Xpsi*B1+[T2'*hatC3';hatC3']*Theta12*[I 0]
ouf_str6='[[B1*Dqw B1*Dqp]+B3]';

% A1'*Xpsi*B1
if size(A1,1)~=0
    str{sc}=['lmiterm([',ks,' 1 3 Xpsi],',ouf_str5,...
        ''',',ouf_str6,');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 2 3 Xpsi],',ouf_str5,...
        ''',',ouf_str6,');'];
    eval(str{sc});
    sc=sc+1;
end

% T2'*hatC3'*Theta12*[I 0]
ouf_str4=['[eye(',nb1,') zeros(',nb1,',',nb2,')]'];

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

% hatC3'*Theta12*[I 0]
str{sc}=['lmiterm([',ks,' 2 3 0],hatC3''*Theta12*',ouf_str4,');'];
eval(str{sc});
sc=sc+1;

%% [T1'*hatA*T2 T1'*hatA]
% T1'*hatA*T2
if size(hatA11,1)~=0
    str{sc}=['lmiterm([',ks,' 4 1 T11],',ouf_str1,...
        '''*hatA11,',ouf_str1,');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 4 1 0],',ouf_str1,...
        '''*hatA12*',ouf_str2,');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 4 1 T12],-',ouf_str1,...
        '''*hatA11,',ouf_str2,');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 4 1 T12],',ouf_str1,...
        ''',hatA22*',ouf_str2,');'];
    eval(str{sc});
    sc=sc+1;
    
end

str{sc}=['lmiterm([',ks,' 4 1 T22],',ouf_str2,...
    ''',hatA22*',ouf_str2,');'];
eval(str{sc});
sc=sc+1;

% T1'*hatA
if size(hatA11,1)~=0
    str{sc}=['lmiterm([',ks,' 4 2 0],',ouf_str1,...
        '''*',ouf_str1,'*hatA);'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 4 2 T12],',ouf_str1,...
        ''',',ouf_str2,'*hatA);'];
    eval(str{sc});
    sc=sc+1;
end
str{sc}=['lmiterm([',ks,' 4 2 T22],',ouf_str2,...
    ''',',ouf_str2,'*hatA);'];
eval(str{sc});
sc=sc+1;

%% T1'*[hatB1 hatB2]
if size(hatA11,1)~=0
    str{sc}=['lmiterm([',ks,' 4 3 0],',ouf_str1,...
        '''*',ouf_str1,'*[hatB1 hatB2]);'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 4 3 T12],',ouf_str1,...
        ''',',ouf_str2,'*[hatB1 hatB2]);'];
    eval(str{sc});
    sc=sc+1;
end
str{sc}=['lmiterm([',ks,' 4 3 T22],',ouf_str2,...
    ''',',ouf_str2,'*[hatB1 hatB2]);'];
eval(str{sc});
sc=sc+1;

%% B1'*Xpsi*B1+[D3' [I;0]]*Theat0*[D3;[I 0]]

% B1'*Xpsi*B1
if size(A1,1)~=0
    str{sc}=['lmiterm([',ks,' 3 3 Xpsi],',ouf_str6,...
        ''',',ouf_str6,');'];
    eval(str{sc});
    sc=sc+1;
end

% [I;0]*Theta12'*D3
str{sc}=['lmiterm([',ks,' 3 3 0],(',ouf_str4,...
    '''*Theta12''*hatD3)+(',ouf_str4,...
    '''*Theta12''*hatD3)'');'];
eval(str{sc});
sc=sc+1;

% [I;0]*Theta22*[I 0]
str{sc}=['lmiterm([',ks,' 3 3 min_gamma],-',ouf_str4,''',',ouf_str4,');'];
eval(str{sc});
sc=sc+1;

%% Theta11^(1/2)*[hatC3*T2 hatC3]
if size(hatA11,1)~=0
    str{sc}=['lmiterm([',ks,' 5 1 T11],sqrtm(Theta11)*hatC3*',ouf_str1,...
        ''',',ouf_str1,');'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['lmiterm([',ks,' 5 1 T12],-sqrtm(Theta11)*hatC3*',ouf_str1,...
        ''',',ouf_str2,');'];
    eval(str{sc});
    sc=sc+1;
end
str{sc}=['lmiterm([',ks,' 5 1 0],sqrtm(Theta11)*hatC3*',ouf_str2,...
    '''*',ouf_str2,');'];
eval(str{sc});
sc=sc+1;

str{sc}=['lmiterm([',ks,' 5 2 0],sqrtm(Theta11)*hatC3);'];
eval(str{sc});
sc=sc+1;

%% Theta11^(1/2)*hatD3
str{sc}=['lmiterm([',ks,' 5 3 0],sqrtm(Theta11)*hatD3);'];
eval(str{sc});
sc=sc+1;

%% -T3
if size(hatA11,1)~=0
    str{sc}=['lmiterm([',ks,' 4 4 T11],-',ouf_str1,...
        ''',',ouf_str1,');'];
    eval(str{sc});
    sc=sc+1;
end
str{sc}=['lmiterm([',ks,' 4 4 T22],-',ouf_str2,...
    ''',',ouf_str2,');'];
eval(str{sc});
sc=sc+1;


%% -I
str{sc}=['lmiterm([',ks,' 5 5 0],-eye(',nc3,'));'];
eval(str{sc});
sc=sc+1;

%% Mpi
for j=1:E.nsimple,
    js=num2str(j);
    nh=0;
    ng=0;
    n=E.simples(j);
    for i=1:size(E.X{n},2),
        is=num2str(i);
        vrs = num2str(E.X{n}(1,i));
        
        str{sc}=['lmiterm([',ks,' 1 1 ',vrs,...
            '],hatC12''*MpiL_',js,'_',is,...
            ',MpiR_',js,'_',is,'*hatC12,''s'');'];
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
        
        str{sc}=['lmiterm([',ks,' 3 3 ',vrs,...
            '],hatD12''*MpiL_',js,'_',is,...
            ',MpiR_',js,'_',is,'*hatD12,''s'');'];
        eval(str{sc})
        sc=sc+1;
        
        nh=nh+E.X{n}(3,i);
        ng=ng+E.X{n}(2,i);
    end
end

%%
ny1=num2str(size(hatB',1));
nx1=num2str(size(hatB1,2));
nx2=num2str(size(hatB',2));
nx3=num2str(size(Cv,1));

str{sc}='NWl=[hatE2''*Theta12 hatB'' hatE2''*sqrtm(Theta11)];';
eval(str{sc});
sc=sc+1;

str{sc}='Wl=null(NWl,''r'');';
eval(str{sc});
sc=sc+1;

str{sc}=['Wl1=Wl(1:',nx1,',:);'];
eval(str{sc});
sc=sc+1;

str{sc}=['Wl2=Wl(',cal_str(nx1,'1'),':',cal_str(nx1,nx2),',:);'];
eval(str{sc});
sc=sc+1;

str{sc}=['Wl3=Wl(',cal_str(nx1,nx2,'1'),':end,:);'];
eval(str{sc});
sc=sc+1;

nx1=num2str(size(Wl1,2));
nx2=num2str(2*size(hatA,2));
nx3=num2str(size(hatB2,2));
ny1=nx2;
ny2=num2str(size(Wl1,1));
ny3=nx3;
ny4=num2str(size(Wl2,1));
ny5=num2str(size(Wl3,1));

str{sc}=['ouf_main_left=[zeros(',ny1,',',nx1,') eye(',nx2,...
    ') zeros(',ny1,',',nx3,');',...
    'Wl1 zeros(',ny2,',',cal_str(nx2,nx3),');',...
    'zeros(',ny3,',',cal_str(nx1,nx2),') eye(',nx3,');',...
    'Wl2 zeros(',ny4,',',cal_str(nx2,nx3),');',...
    'Wl3 zeros(',ny5,',',cal_str(nx2,nx3),')];'];
eval(str{sc});
sc=sc+1;

str{sc}=['lmiterm([',ks,' 0 0 0],ouf_main_left);'];
eval(str{sc})
sc=sc+1;

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

if ~isinf(gain)
    
    str{sc}='setlmis([]);';
    eval(str{sc});
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
    
    str{sc}=['T2=[T11 -T12;zeros(',nx2,',',nx1,') eye(',nx2,')];'];
    eval(str{sc})
    sc=sc+1;
    
    if size(A1,1)~=0
        str{sc}='Xpsi=decs2mats(xopt,sXpsi);';
        eval(str{sc})
        sc=sc+1;
        
        str{sc}=['Xpsie=',ouf_str3,'''*decs2mats(xopt,sXpsi)*',ouf_str3,';'];
        eval(str{sc})
        sc=sc+1;
    else
        str{sc}='Xpsi=[];';
        eval(str{sc})
        sc=sc+1;
        
        str{sc}=['Xpsie=zeros(',num2str(size(X,1)),');'];
        eval(str{sc})
        sc=sc+1;
    end
    
    % define M_E
    nx1=num2str(size(hatA,1)+size(Cy,1));
    ny1=num2str(size(hatA,1)+size(Cv,1));
    str{sc}=['[M_E,nM_E,sM_E]=lmivar(2,[',ny1,' ',nx1,']);'];
    eval(str{sc});
    sc=sc+1;
    
    nx1=num2str(size(hatA11,1));
    nx2=num2str(size(hatA22,1));
    
    str{sc}=['bY=[[T11 zeros(',nx1,',',nx2,...
        ');zeros(',nx2,',',nx1,') T22] T1'';T1 X];'];
    eval(str{sc});
    sc=sc+1;
    
    nx1=num2str(size(T2,2));
    nx2=num2str(size(hatA,2));
    nx3=num2str(size(hatB1,2));
    nx4=num2str(size(hatB2,2));
    nx5=num2str(size(bY,2));
    
    ny2=num2str(size(A1,1));
    ny3=num2str(size(T1',1));
    ny4=num2str(size(X,1));
    ny5=num2str(size(hatC12,1));
    ny6=num2str(size(hatC3,1));
    
    str{sc}=['ouf_Xi11=[eye(',cal_str(nx1,nx2),') zeros(',cal_str(nx1,nx2),...
        ',',cal_str(nx3,nx4,nx5),');',...
        ouf_str5,' ',ouf_str5,' ',ouf_str6,' zeros(',ny2,',',nx5,');',...
        'T1''*hatA*T2 T1''*hatA T1''*hatB1 T1''*hatB2 zeros(',ny3,',',nx5,');',...
        'zeros(',ny4,',',nx1,') X*hatA X*hatB1 X*hatB2 zeros(',ny4,',',nx5,');',...
        'zeros(',nx5,',',cal_str(nx1,nx2,nx3,nx4),') eye(',nx5,');',...
        'hatC12 hatC12 hatD12 zeros(',ny5,',',nx5,');',...
        'hatC3*T2 hatC3 hatD3 zeros(',ny6,',',nx5,');',...
        'zeros(',nx3,',',cal_str(nx1,nx2),') eye(',nx3,') zeros(',nx3,...
        ',',cal_str(nx4,nx5),')];'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['Xi21=sqrtm(Theta11)*[hatC3*T2 hatC3 hatD3 zeros(',...
        ny6,',',nx5,')];'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}=['Xi22=-eye(',ny6,');'];
    eval(str{sc});
    sc=sc+1;
    
    nx6=ny6;
    ny1=num2str(size(hatA,1));
    ny2=num2str(size(Cy,1));
    str{sc}=['mathC=[eye(',ny1,') zeros(',ny1,',',...
        cal_str(nx2,nx3,nx4,nx5,nx6),');',...
        'zeros(',ny2,',',nx1,') hatC hatF1 hatF2 zeros(',ny2,...
        ',',cal_str(nx5,nx6),')];'];
    eval(str{sc});
    sc=sc+1;
    
    ny1=num2str(size(hatA,1));
    ny2=num2str(size(Cv,1));
    nx1=num2str(2*size(hatA,2));
    nx2=num2str(size(hatB1,2));
    nx3=num2str(size(hatB2,2));
    nx4=num2str(size(hatB',2));
    nx5=nx4;
    nx6=num2str(size(Cv,1));
    str{sc}=['mathB=[zeros(',ny1,',',cal_str(nx1,nx2,nx3,nx4),') eye(',nx5,...
        ') zeros(',ny1,',',nx6,');',...
        'zeros(',ny2,',',nx1,') hatE2''*Theta12 zeros(',ny2,...
        ',',nx3,') hatB'' zeros(',ny2,',',...
        nx5,') hatE2''*sqrtm(Theta11)];'];
    eval(str{sc});
    sc=sc+1;
    
    % Mpi
    nx1=num2str(size(hatC12,1));
    str{sc}=['Mpi=zeros(',nx1,',',nx1,');'];
    eval(str{sc});
    sc=sc+1;
    for j=1:E.nsimple,
        js=num2str(j);
        nh=0;
        ng=0;
        n=E.simples(j);
        for i=1:size(E.X{n},2),
            is=num2str(i);
            vr = E.X{n}(1,i);
            if vr>0
                str{sc}=['Mpi=Mpi+MpiL_',js,'_',is,'*decs2mats(xopt,sX',...
                    num2str(vr),')*MpiR_',js,'_',is,...
                    '+(MpiL_',js,'_',is,'*decs2mats(xopt,sX',...
                    num2str(vr),')*MpiR_',js,'_',is,')'';']; %#ok<*AGROW>
            else
                str{sc}=['Mpi=Mpi+MpiL_',js,'_',is,'*decs2mats(xopt,sX',...
                    num2str(-vr),')''*MpiR_',js,'_',is,...
                    '+(MpiL_',js,'_',is,'*decs2mats(xopt,sX',...
                    num2str(-vr),')''*MpiR_',js,'_',is,')'';'];
            end
            eval(str{sc});
            sc=sc+1;
            nh=nh+E.X{n}(3,i);
            ng=ng+E.X{n}(2,i);
        end
    end
    
    nx1=num2str(size(Mpi,1));
    nx2=num2str(size(Theta11,1));
    nx3=num2str(size(Theta12,2));
    str{sc}=['Pr0=[Mpi zeros(',nx1,',',cal_str(nx2,nx3),');',...
        'zeros(',nx2,',',cal_str(nx1,nx2),') Theta12;',...
        'zeros(',nx3,',',nx1,') Theta12'' -',num2str(gain^2),...
        '*eye(',nx3,')];'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='Xclpsie=[Xpsie Xpsie;Xpsie Xpsie];';
    eval(str{sc});
    sc=sc+1;
    
    nx1=num2str(size(bY,2));
    nx2=num2str(size(Xpsi,2));
    nx3=nx1;
    nx4=nx1;
    nx5=num2str(size(Pr0,2));
    str{sc}=['M_Xi11=[-(bY+Xclpsie) zeros(',nx1,',',...
        cal_str(nx2,nx3,nx4,nx5),');',...
        'zeros(',nx2,',',nx1,') Xpsi zeros(',nx2,',',...
        cal_str(nx3,nx4,nx5),');',...
        'zeros(',nx3,',',cal_str(nx1,nx2,nx3),') eye(',nx4,...
        ') zeros(',nx3,',',nx5,');',...
        'zeros(',nx4,',',cal_str(nx1,nx2),') eye(',nx3,...
        ') -bY zeros(',nx4,',',nx5,');',...
        'zeros(',nx5,',',cal_str(nx1,nx2,nx3,nx4),') Pr0];'];
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='Xi11=ouf_Xi11''*M_Xi11*ouf_Xi11;';
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='Xi11=(Xi11+Xi11'')/2;';
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='Xi=[Xi11 Xi21'';Xi21 Xi22];';
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='lmiterm([1 1 1 0],Xi);';
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='lmiterm([1 1 1 M_E],mathB'',mathC,''s'');';
    eval(str{sc});
    sc=sc+1;
    
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
    
    str{sc}='[feas,xopt]=feasp(lmi,lmiparameter);';
    eval(str{sc});
    sc=sc+1;
    
    if feas<0
        if vrb,
            disp_str(67,'feasible');
        end
        feas=1;
    else
        if vrb,
            disp_str(67,'infeasible');
        end
        feas=0;
    end
    
    estimation=[];
    
    if feas
        
        str{sc}='M_E=decs2mats(xopt,sM_E);';
        eval(str{sc});
        sc=sc+1;
        
        ny1=num2str(size(hatA,1));
        
        str{sc}=['hatK=M_E(1:',ny1,',1:',ny1,');'];
        eval(str{sc});
        sc=sc+1;
        
        str{sc}=['L=M_E(1:',ny1,',',ny1,'+1:end);'];
        eval(str{sc});
        sc=sc+1;
        
        str{sc}=['hatM=M_E(',ny1,'+1:end,1:',ny1,');'];
        eval(str{sc});
        sc=sc+1;
        
        str{sc}=['N=M_E(',ny1,'+1:end,',ny1,'+1:end);'];
        eval(str{sc});
        sc=sc+1;
        
        getEST;
        if any(abs(eig(AE))>1)
            disp_str(59,'estimation unstable')
            return
        end
        
        str{sc}='estimation=ss(AE,BE,CE,DE,-1);';
        eval(str{sc});
        sc=sc+1;
    end
end

if lmiparameter(5)>100;
    save iqc_estimation_l2gain_DTeli_lmilab_exe A A1 A2 A3 B1 B2 B3...
        Dqp Dqw Bp Bw C1 C2 D1 D2 Cq Cv Cy Dvp Dvw Dyp Dyw...
        lmiparameter A_W B_W C_W D_W
    if exist('state1','var')
        save iqc_estimation_l2gain_DTeli_lmilab_exe state* -append
    end
    if exist('struct1','var')
        save iqc_estimation_l2gain_DTeli_lmilab_exe struct* -append
    end
    if exist('l_1_1','var')
        save iqc_estimation_l2gain_DTeli_lmilab_exe l_* -append
    end
    if exist('r_1_1','var')
        save iqc_estimation_l2gain_DTeli_lmilab_exe r_* -append
    end
    if exist('MpiL_1_1','var')
        save iqc_estimation_l2gain_DTeli_lmilab_exe MpiL_* -append
    end
    if exist('MpiR_1_1','var')
        save iqc_estimation_l2gain_DTeli_lmilab_exe MpiR_* -append
    end
    fid=fopen('iqc_estimation_l2gain_DTeli_lmilab_exe.m','wt');
    for i=1:sc-1
        fprintf(fid,[str{i} '\n']);
    end
    fclose(fid);
end
