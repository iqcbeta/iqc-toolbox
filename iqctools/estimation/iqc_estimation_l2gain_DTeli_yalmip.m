function [gain,estimation]=...
    iqc_estimation_l2gain_DTeli_yalmip(w,v,y,p,q,weight)
% function [gain,estimation]=...
%     iqc_estimation_l2gain_DTeli_yalmip(w,v,y,p,q,weight)
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
sft_v='(1e-8)';

switch ABST.systemtype
    case 'continuous',
        disp_str(70,'iqc_estimation_l2gain_DTeli_yalmip','continuous')
    case 'discrete',
        psi='[1 0;0 -1]';
end
num_var_pos={};

lmiparameter=ABST.lmiparameter;

% vrb=0 to suppress the output
vrb=1;
if ~isempty(ABST.lmiparameter)
    vrb = isempty(strfind(ABST.lmiparameter,'''verbose'',0'));
end

if ~vrb
    disp_str(50,'iqc_estimation_l2gain_DTeli_yalmip')
else
    disp_str(51,'iqc_estimation_l2gain_DTeli_yalmip')
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
str{1}=['%% Created by iqc_estimation_l2gain_DTeli_yalmip on ' export_date];
str{2}='\n%% Load coefficient matrices ...';
str{3}='load iqc_estimation_l2gain_DTeli_yalmip_exe';

ALL_LMI=[];
%% initialize the YALMIP LMI
str{4}='\n%% Initialize the YALMIP LMI';
str{5}='ALL_LMI = [];';
eval(str{5});

%% variable definite
defVar;

%% now we define the "non-KYP" lmi
nonKYP;

%% MpiR_*_* MpiL_*_*
multiplier_M;

%% Mpi
Mpi=zeros(size(MpiR_1_1,2));
for j=1:E.nsimple,
    js=num2str(j);
    nh=0;
    ng=0;
    n=E.simples(j);
    for i=1:size(E.X{n},2),
        is=num2str(i);
        vr = E.X{n}(1,i);
        
        if vr < 0;
            vrs=num2str(-vr);
            str{sc}=['Mpi=Mpi+MpiL_',js,'_',is,'*X',vrs,...
                '''*MpiR_',js,'_',is,'+(MpiL_',js,'_',is,'*X',vrs,...
                '''*MpiR_',js,'_',is,')'';']; %#ok<*AGROW>
            eval(str{sc});
            sc=sc+1;
        else
            vrs=num2str(vr);
            str{sc}=['Mpi=Mpi+MpiL_',js,'_',is,'*X',vrs,...
                '*MpiR_',js,'_',is,'+(MpiL_',js,'_',is,'*X',vrs,...
                '*MpiR_',js,'_',is,')'';'];
            eval(str{sc});
            sc=sc+1;
        end
        nh=nh+E.X{n}(3,i);
        ng=ng+E.X{n}(2,i);
    end
end

%% Mpi11>0
Mpi11gt0;

%% now we define the "main" lmi
if vrb
    disp_str(54)
end

str{sc}='\n%% define the "main" lmi terms...';
sc=sc+1;

str{sc}='min_gamma=sdpvar;';
eval(str{sc});
sc=sc+1;
num_var_pos=[num_var_pos;{'min_gamma','1',[]}];

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

str{sc}=['X=sdpvar(',num2str(size(hatA,1)),');'];
eval(str{sc});
sc=sc+1;
num_var_pos=[num_var_pos;{'X','1',[]}];

nx1=num2str(size(hatA11,1));
nx2=num2str(size(hatA22,1));
str{sc}=['T11=sdpvar(',nx1,');'];
eval(str{sc});
sc=sc+1;
num_var_pos=[num_var_pos;{'T11','1',[]}];

str{sc}=['T12=sdpvar(',nx1,',',nx2,',''full'');'];
eval(str{sc});
sc=sc+1;
num_var_pos=[num_var_pos;{'T12','2',[]}];

str{sc}=['T22=sdpvar(',nx2,');'];
eval(str{sc});
sc=sc+1;
num_var_pos=[num_var_pos;{'T22','1',[]}];

str{sc}=['T1=[eye(',nx1,') zeros(',nx1,',',nx2,');',...
    'T12'' T22];'];
eval(str{sc});
sc=sc+1;

str{sc}=['T2=[T11 -T12;zeros(',nx2,',',nx1,') eye(',nx2,')];'];
eval(str{sc});
sc=sc+1;

str{sc}=['T3=[T11 zeros(',nx1,',',nx2,');zeros(',nx2,',',nx1,') T22];'];
eval(str{sc});
sc=sc+1;

str{sc}='Ycl=[T3 T1'';T1 X];';
eval(str{sc});
sc=sc+1;

nx1=num2str(size(Ycl,1));
str{sc}=['ALL_LMI=ALL_LMI+[Ycl>=',sft_v,'*eye(',nx1,')];'];
eval(str{sc});
sc=sc+1;

str{sc}=['Acl1=[zeros(',na2,',',na1,') A1 A3 B1*Cq',...
    ' zeros(',na2,',',na1,') A1 A3 B1*Cq];'];
eval(str{sc});
sc=sc+1;

str{sc}='Bcl1=[B1*Dqw B1*Dqp]+B3;';
eval(str{sc});
sc=sc+1;

nx1=na1;
nx2=num2str(size(T22,1));
str{sc}=['T1AT2=[hatA11*T11 hatA12-hatA11*T12+T12*hatA22;',...
    'zeros(',nx2,',',nx1,') T22*hatA22];'];
eval(str{sc});
sc=sc+1;

str{sc}=['Theta11=eye(',nc3,');'];
eval(str{sc});
sc=sc+1;

str{sc}=['Theta12=zeros(',nc3,',',nb1,');'];
eval(str{sc});
sc=sc+1;

nx1=num2str(size(Mpi,1));
nx2=num2str(size(Cv,1));
nx3=num2str(size(Bw,2));
str{sc}=['Pr=[Mpi zeros(',nx1,',',cal_str(nx2,nx3),');',...
    'zeros(',nx2,',',nx1,') Theta11 Theta12;',...
    'zeros(',nx3,',',nx1,') Theta12'' -min_gamma*eye(',nx3,')];'];
eval(str{sc});
sc=sc+1;
str{sc}=['Pr0=[Mpi zeros(',nx1,',',cal_str(nx2,nx3),');',...
    'zeros(',nx2,',',cal_str(nx1,nx2),') Theta12;',...
    'zeros(',nx3,',',nx1,') Theta12'' -min_gamma*eye(',nx3,')];'];
eval(str{sc});
sc=sc+1;

str{sc}=['Xpsie=[zeros(',na1,',',cal_str(na1,na2,na3,na4),');',...
    'zeros(',na2,',',na1,') Xpsi zeros(',na2,',',cal_str(na3,na4),');',...
    'zeros(',cal_str(na3,na4),',',cal_str(na1,na2,na3,na4),')];'];
eval(str{sc});
sc=sc+1;

str{sc}='Xclpsie=[Xpsie Xpsie;Xpsie Xpsie];';
eval(str{sc});
sc=sc+1;

str{sc}='\n%% define the "main" right lmi terms...';
sc=sc+1;

nx1=num2str(size(hatA,2));
nx2=num2str(size(hatB1,2));
nx3=num2str(size(hatB2,2));
str{sc}=['ouf_main_right=[eye(',nx1,') zeros(',nx1,',',cal_str(nx2,nx3),');',...
    'hatA hatB1 hatB2;',...
    'hatC12 hatD12;',...
    'hatC3 hatD3;',...
    'zeros(',nx2,',',nx1,') eye(',nx2,') zeros(',nx2,',',nx3,')];'];
eval(str{sc});
sc=sc+1;

nx1=num2str(size(X,2));
nx2=nx1;
nx3=num2str(size(Pr,2));
str{sc}=['M_main_right=[-(X+Xpsie) zeros(',nx1,',',cal_str(nx2,nx3),');',...
    'zeros(',nx2,',',nx1,') (X+Xpsie) zeros(',nx2,',',nx3,');',...
    'zeros(',nx3,',',cal_str(nx1,nx2),') Pr];'];
eval(str{sc});
sc=sc+1;

str{sc}='NWr=[hatC hatF1 hatF2];';
eval(str{sc});
sc=sc+1;

str{sc}='Wr=null(NWr,''r'');';
eval(str{sc});
sc=sc+1;

nx1=num2str(size(Wr,2));
str{sc}=['ALL_LMI=ALL_LMI+',...
    '[Wr''*ouf_main_right''*M_main_right*ouf_main_right*Wr',...
    '<=-',sft_v,'*eye(',nx1,')];'];
eval(str{sc});
sc=sc+1;

str{sc}='\n%% define the "main" left lmi terms...';
sc=sc+1;

nx1=num2str(size(Ycl,2));
nx2=num2str(size(Xpsi,2));
nx3=num2str(size(T3,2));
nx4=nx3;
nx5=num2str(size(Pr0,2));
str{sc}=['M_Gamma_left11=[-(Ycl+Xclpsie) zeros(',nx1,',',...
    cal_str(nx2,nx3,nx4,nx5),');',...
    'zeros(',nx2,',',nx1,') Xpsi zeros(',nx2,',',...
    cal_str(nx3,nx4,nx5),');',...
    'zeros(',nx3,',',cal_str(nx1,nx2,nx3),') eye(',nx4,...
    ') zeros(',nx3,',',nx5,');',...
    'zeros(',nx4,',',cal_str(nx1,nx2),') eye(',nx3,...
    ') -T3 zeros(',nx4,',',nx5,');',...
    'zeros(',nx5,',',cal_str(nx1,nx2,nx3,nx4),') Pr0];'];
eval(str{sc});
sc=sc+1;

nx1=num2str(size(hatA,2));
nx2=nx1;
nx3=num2str(size(hatB1,2));
nx4=num2str(size(hatB2,2));
nx5=num2str(size(T3,2));

ny3=num2str(size(Acl1,1));
ny4=num2str(size(T3,1));
ny6=num2str(size(hatC12,1));
ny7=num2str(size(hatC3,1));

str{sc}=['ouf_Gamma_left11=[eye(',cal_str(nx1,nx2),...
    ') zeros(',cal_str(nx1,nx2),',',cal_str(nx3,nx4,nx5),');',...
    'Acl1 Bcl1 zeros(',ny3,',',nx5,');',...
    'T1AT2 T1''*hatA T1''*hatB1 T1''*hatB2 zeros(',ny4,',',nx5,');',...
    'zeros(',nx5,',',cal_str(nx1,nx2,nx3,nx4),') eye(',nx5,');',...
    'hatC12 hatC12 hatD12 zeros(',ny6,',',nx5,');',...
    'hatC3*T2 hatC3 hatD3 zeros(',ny7,',',nx5,');',...
    'zeros(',nx3,',',cal_str(nx1,nx2),') eye(',nx3,') zeros(',nx3,...
    ',',cal_str(nx4,nx5),')];'];
eval(str{sc});
sc=sc+1;

str{sc}='Gamma_left11=ouf_Gamma_left11''*M_Gamma_left11*ouf_Gamma_left11;';
eval(str{sc});
sc=sc+1;

nx1=num2str(size(T3,1));
ny1=num2str(size(hatC3,1));
str{sc}=['Gamma_left21=sqrtm(Theta11)*[hatC3*T2 hatC3 hatD3 zeros(',...
    ny1,',',nx1,')];'];
eval(str{sc});
sc=sc+1;

str{sc}=['Gamma_left=[Gamma_left11 Gamma_left21'';',...
    'Gamma_left21 -eye(',ny1,')];'];
eval(str{sc});
sc=sc+1;

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

nx1=num2str(size(ouf_main_left,2));
str{sc}=['ALL_LMI=ALL_LMI+[ouf_main_left''*Gamma_left*ouf_main_left',...
    '<=-',sft_v,'*eye(',nx1,')];'];
eval(str{sc});
sc=sc+1;

%% get variable
xnumber=[];
for i1=1:size(num_var_pos,1)
    eval(['[di_y,di_x]=size(',num_var_pos{i1,1},');']);
    switch num_var_pos{i1,2}
        case '1'
            for i2=1:di_x
                eval(['newxnumber=',num_var_pos{i1,1},'(',...
                    '1:',num2str(i2),',',num2str(i2),');']);
                newxnumber=reshape(newxnumber,1,size(newxnumber,1));
                xnumber=[xnumber,newxnumber];
            end
        case '2'
            for i2=1:di_y
                eval(['newxnumber=',num_var_pos{i1,1},'(',num2str(i2),...
                    ',:);']);
                xnumber=[xnumber,newxnumber];
            end
        case '3';
            for i2=1:di_x
                eval(['xnumber=[xnumber,',num_var_pos{i1,1},'(',...
                    num2str(i2),',',num2str(i2),')];']);
            end
        case '4';
            for i2=1:(di_x-1)
                eval(['xnumber=[xnumber,',num_var_pos{i1,1},'(',...
                    num2str(i2),',',num2str(i2+1),':',...
                    num2str(di_x),')];']);
            end
        case '5';
            for i2=1:size(num_var_pos{i1,3},1)
                eval(['newxnumber=num_var_pos{i1,3}(',...
                    num2str(i2),',:);']);
                eval(['xnumber=[xnumber,',newxnumber,'];']);
            end
    end
end
xnumber=reshape(xnumber,size(xnumber,2),1);

%% solving the system of LMI's
str{sc}='\n%% Solving the system of LMIs ...';
sc=sc+1;
if ~isempty(ABST.lmiparameter),
    str{sc}=['ops=sdpsettings(',ABST.lmiparameter,');'];
else
    str{sc}='ops=[];';
end
eval(str{sc});
sc=sc+1;

if vrb,
    disp_str(56,num2str(length(xnumber)))
end

str{sc}='sol = solvesdp(ALL_LMI,min_gamma,ops);';
eval(str{sc});
sc=sc+1;

gain=inf;

checkset(ALL_LMI)
switch sol.problem
    case 0 % success
        if vrb
            disp_str(57,'L2-gain','Feasible')
        end
        xopt=double(xnumber);
        str{sc} = 'gain=sqrt(double(min_gamma))';
        eval([str{sc},';'])
        sc=sc+1;
    case 1
        if vrb
            disp_str(57,'L2-gain','Infeasible')
        end
        xopt=[];
        str{sc} = 'gain=inf';
        eval([str{sc},';'])
        sc=sc+1;
    case {4,5}
        alleig=[];
        for i1=1:size(ALL_LMI)
            alleig=[alleig;eig(double(ALL_LMI{i1}+...
                eval(sft_v)*eye(size(ALL_LMI{i1},1))))];
        end
        if min(alleig)>0
            if vrb
                disp_str(57,'L2-gain','Feasible');
            end
            xopt=double(xnumber);
            str{sc} = 'gain=sqrt(double(min_gamma))';
            eval([str{sc},';'])
            sc=sc+1;
        elseif min(alleig)<0 && -min(alleig) < str2double(sft_v)
            if vrb
                disp_str(58,num2str(min(alleig)))
            end
            xopt=double(xnumber);
            str{sc} = 'gain=sqrt(double(min_gamma))';
            eval([str{sc},';'])
            sc=sc+1;
        else
            if vrb
                disp_str(57,'L2-gain','Infeasible')
            end
            xopt=[];
            str{sc} = 'gain=inf';
            eval([str{sc},';'])
            sc=sc+1;
        end
    otherwise
        if vrb
            disp_str(59,yalmiperror(sol.problem))
        end
        xopt=[];
        str{sc} = 'gain=inf';
        eval([str{sc},';'])
        sc=sc+1;
end

%% ­««Ø¦ô´ú¾¹
estimation=[];

if ~isinf(gain)
    
    num_var_pos={};
    
    str{sc}='ALL_LMI = [];';
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='X=double(X);';
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='T11=double(T11);';
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='T12=double(T12);';
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='T22=double(T22);';
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='T1=double(T1);';
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='T2=double(T2);';
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='Ycl=double(Ycl);';
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='Pr0=double(Pr0);';
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='Xpsi=double(Xpsi);';
    eval(str{sc});
    sc=sc+1;
    
    str{sc}='Xclpsie=double(Xclpsie);';
    eval(str{sc});
    sc=sc+1;
    
    % define hatM hatK L N
    nx1=num2str(size(hatA,1));
    str{sc}=['hatK=sdpvar(',nx1,',',nx1,',''full'');'];
    eval(str{sc});
    sc=sc+1;
    num_var_pos=[num_var_pos;{'hatK','2',[]}];
    
    nx2=num2str(size(Cy,1));
    str{sc}=['L=sdpvar(',nx1,',',nx2,',''full'');'];
    eval(str{sc});
    sc=sc+1;
    num_var_pos=[num_var_pos;{'L','2',[]}];
    
    ny2=num2str(size(Cv,1));
    str{sc}=['hatM=sdpvar(',ny2,',',nx1,',''full'');'];
    eval(str{sc});
    sc=sc+1;
    num_var_pos=[num_var_pos;{'hatM','2',[]}];
    
    str{sc}=['N=sdpvar(',ny2,',',nx2,',''full'');'];
    eval(str{sc});
    sc=sc+1;
    num_var_pos=[num_var_pos;{'N','2',[]}];
    
    str{sc}='M_E=[hatK L;hatM N];';
    eval(str{sc});
    sc=sc+1;
    
    nx1=num2str(size(T2,2));
    nx2=num2str(size(hatA,2));
    nx3=num2str(size(hatB1,2));
    nx4=num2str(size(hatB2,2));
    nx5=num2str(size(Ycl,2));
    
    ny2=num2str(size(Acl1,1));
    ny3=num2str(size(T1',1));
    ny4=num2str(size(X,1));
    ny5=num2str(size(hatC12,1));
    ny6=num2str(size(hatC3,1));
    
    str{sc}=['ouf_Xi11=[eye(',cal_str(nx1,nx2),') zeros(',cal_str(nx1,nx2),...
        ',',cal_str(nx3,nx4,nx5),');',...
        'Acl1 Bcl1 zeros(',ny2,',',nx5,');',...
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
    
    nx1=num2str(size(Ycl,2));
    nx2=num2str(size(Xpsi,2));
    nx3=nx1;
    nx4=nx1;
    nx5=num2str(size(Pr0,2));
    str{sc}=['M_Xi11=[-(Ycl+Xclpsie) zeros(',nx1,',',...
        cal_str(nx2,nx3,nx4,nx5),');',...
        'zeros(',nx2,',',nx1,') Xpsi zeros(',nx2,',',...
        cal_str(nx3,nx4,nx5),');',...
        'zeros(',nx3,',',cal_str(nx1,nx2,nx3),') eye(',nx4,...
        ') zeros(',nx3,',',nx5,');',...
        'zeros(',nx4,',',cal_str(nx1,nx2),') eye(',nx3,...
        ') -Ycl zeros(',nx4,',',nx5,');',...
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
    
    str{sc}='F=Xi+(mathB''*M_E*mathC)+(mathB''*M_E*mathC)'';';
    eval(str{sc});
    sc=sc+1;
    
    nx1=num2str(size(F,2));
    str{sc}=['ALL_LMI=ALL_LMI+[F<=-',sft_v,'*eye(',nx1,')];'];
    eval(str{sc});
    sc=sc+1;
    
    %% get variable
    xnumber=[];
    for i1=1:size(num_var_pos,1)
        eval(['[di_y,di_x]=size(',num_var_pos{i1,1},');']);
        switch num_var_pos{i1,2}
            case '1'
                for i2=1:di_x
                    eval(['newxnumber=',num_var_pos{i1,1},'(',...
                        '1:',num2str(i2),',',num2str(i2),');']);
                    newxnumber=reshape(newxnumber,1,size(newxnumber,1));
                    xnumber=[xnumber,newxnumber];
                end
            case '2'
                for i2=1:di_y
                    eval(['newxnumber=',num_var_pos{i1,1},'(',num2str(i2),...
                        ',:);']);
                    xnumber=[xnumber,newxnumber];
                end
            case '3';
                for i2=1:di_x
                    eval(['xnumber=[xnumber,',num_var_pos{i1,1},'(',...
                        num2str(i2),',',num2str(i2),')];']);
                end
            case '4';
                for i2=1:(di_x-1)
                    eval(['xnumber=[xnumber,',num_var_pos{i1,1},'(',...
                        num2str(i2),',',num2str(i2+1),':',...
                        num2str(di_x),')];']);
                end
            case '5';
                for i2=1:size(num_var_pos{i1,3},1)
                    eval(['newxnumber=num_var_pos{i1,3}(',...
                        num2str(i2),',:);']);
                    eval(['xnumber=[xnumber,',newxnumber,'];']);
                end
        end
    end
    xnumber=reshape(xnumber,size(xnumber,2),1);
    
    %% solving the system of LMI's
    str{sc}='\n%% Solving the system of LMIs ...';
    sc=sc+1;
    
    if vrb,
        disp_str(56,num2str(length(xnumber)))
    end
    
    str{sc}='sol = solvesdp(ALL_LMI,[],ops);';
    eval(str{sc});
    sc=sc+1;
    
    checkset(ALL_LMI)
    switch sol.problem
        case 0 % success
            if vrb
                disp_str(67,'feasible')
            end
            feas=1;
        case 1
            if vrb
                disp_str(67,'infeasible')
            end
            feas=0;
        case {4,5}
            alleig=[];
            for i1=1:size(ALL_LMI)
                alleig=[alleig;eig(double(ALL_LMI{i1}+...
                    eval(sft_v)*eye(size(ALL_LMI{i1},1))))];
            end
            if min(alleig)>0
                if vrb
                    disp_str(67,'feasible')
                end
                feas=1;
            elseif min(alleig)<0 && -min(alleig) < str2double(sft_v)
                if vrb
                    disp_str(58,num2str(min(alleig)))
                end
                feas=1;
            else
                if vrb
                    disp_str(67,'infeasible')
                end
                feas=0;
            end
        otherwise
            if vrb
                disp_str(59,yalmiperror(sol.problem))
            end
            feas=0;
    end
    
    if feas==1
        ny1=num2str(size(hatA,1));
        str{sc}=['hatK=double(M_E(1:',ny1,',','1:',ny1,'));'];
        eval(str{sc});
        sc=sc+1;
        
        str{sc}=['L=double(M_E(1:',ny1,',','1+',ny1,':end));'];
        eval(str{sc});
        sc=sc+1;
        
        str{sc}=['hatM=double(M_E(1+',ny1,':end,','1:',ny1,'));'];
        eval(str{sc});
        sc=sc+1;
        
        str{sc}=['N=double(M_E(1+',ny1,':end,','1+',ny1,':end));'];
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
    else
        str{sc}='estimation=[];';
        eval(str{sc});
        sc=sc+1;
    end
    
end

if ~vrb
    save iqc_estimation_l2gain_DTeli_yalmip_exe A A1 A2 A3 B1 B2 B3...
        Dqp Dqw Bp Bw C1 C2 D1 D2 Cq Cv Cy Dvp Dvw Dyp Dyw...
        lmiparameter A_W B_W C_W D_W
    if exist('state1','var')
        save iqc_estimation_l2gain_DTeli_yalmip_exe state* -append
    end
    if exist('struct1','var')
        save iqc_estimation_l2gain_DTeli_yalmip_exe struct* -append
    end
    if exist('l_1_1','var')
        save iqc_estimation_l2gain_DTeli_yalmip_exe l_* -append
    end
    if exist('r_1_1','var')
        save iqc_estimation_l2gain_DTeli_yalmip_exe r_* -append
    end
    if exist('MpiL_1_1','var')
        save iqc_estimation_l2gain_DTeli_yalmip_exe MpiL_* -append
    end
    if exist('MpiR_1_1','var')
        save iqc_estimation_l2gain_DTeli_yalmip_exe MpiR_* -append
    end
    fid=fopen('iqc_estimation_l2gain_DTeli_yalmip_exe.m','wt');
    for i=1:sc-1
        fprintf(fid,[str{i} '\n']);
    end
    fclose(fid);
end
