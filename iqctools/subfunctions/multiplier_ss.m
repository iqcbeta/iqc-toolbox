function [Psi]=multiplier_ss(E,q,p,w)
% q: uncertainty input
% p: uncertainty output
% w: exogenous excitation
% Psi: [Psi_q Psi_p Psi_w]

global ABST
log=ABST.log;

%% 假如 q, p, w 使用者輸入的為訊號組合，如 q=[q1;q2]...
%% q
nq = double(q);
q_sgn_vec=[];
q_sgn_vec=sgn_vec(nq,0,log,q_sgn_vec);
q_sgn_vec=unique(q_sgn_vec);

%% p
np = double(p);
p_sgn_vec=[];
p_sgn_vec=sgn_vec(np,0,log,p_sgn_vec);
p_sgn_vec=unique(p_sgn_vec);

%% w
nw = double(w);
w_sgn_vec=[];
w_sgn_vec=sgn_vec(nw,0,log,w_sgn_vec);
w_sgn_vec=unique(w_sgn_vec);

%%

step_sgn=unique([q_sgn_vec,p_sgn_vec,w_sgn_vec]);

table_id = find(E.T(:,1)==11);

M_state_pos=[];
for i1=1:length(table_id)
    
    get_table=disp_table(table_id(i1),step_sgn);
    
    state_pos=[];
    for i2=1:size(get_table,1)
        if length(E.state_lnk_pos)>=get_table(i2)
            state_pos=[state_pos,E.state_lnk_pos{get_table(i2)}]; %#ok<*AGROW>
        end
    end
    state_pos=sort(state_pos);
    M_state_pos=[M_state_pos,state_pos];
end
M_state_pos=unique(M_state_pos);

N_state_pos=setdiff(1:E.nstates,M_state_pos);
%%
Acl=E.ab(:,1:E.nstates);
Bcl=E.ab(:,E.nstates+1:end);

CD_p=E.C{np};
[nx,ny1]=find(CD_p); %#ok<*ASGLU>
pos_p=unique(ny1');
CD_w=E.C{nw};
[nx,ny2]=find(CD_w);
pos_w=unique(ny2');

CD_q=E.C{nq};
C_q=CD_q(:,N_state_pos);
D_qp=CD_q(:,pos_p);
D_qw=CD_q(:,pos_w);

Api=Acl(M_state_pos,M_state_pos);

BpiqCq=Acl(M_state_pos,N_state_pos);

Bpi_q=BpiqCq/C_q;
if ~isempty(find(Bpi_q*C_q-BpiqCq>eps, 1))
    disp_str(72)
end

BpiqDqw_Bpiw=Bcl(M_state_pos,pos_w-E.nstates);
Bpi_w=BpiqDqw_Bpiw-Bpi_q*D_qw;
if ~isempty(find((Bpi_w+Bpi_q*D_qw)-BpiqDqw_Bpiw>eps, 1))
    disp_str(72)
end

BpiqDqp_Bpip=Bcl(M_state_pos,pos_p-E.nstates);
Bpi_p=BpiqDqp_Bpip-Bpi_q*D_qp;
if ~isempty(find((Bpi_p+Bpi_q*D_qp)-BpiqDqp_Bpip>eps, 1))
    disp_str(72)
end

Bpi=[Bpi_q,Bpi_w,Bpi_p];
na=size(Api,1);
nb=size(Bpi,2);
Psi_ss=ss(Api,Bpi,eye(na),zeros(na,nb));
Psi_ss=[Psi_ss;eye(nb)];

pos_q=1:size(Bpi_q,2);
pos_w=pos_q(end)+1:pos_q(end)+size(Bpi_w,2);
pos_p=pos_w(end)+1:size(Bpi,2);

%% 儲存資料
Psi.Api=Psi_ss.a;
Psi.Bpi_q=Psi_ss.b(:,pos_q);
Psi.Bpi_p=Psi_ss.b(:,pos_p);
Psi.Bpi_w=Psi_ss.b(:,pos_w);
Psi.Cpi=Psi_ss.c;
Psi.Dpi_q=Psi_ss.d(:,pos_q);
Psi.Dpi_p=Psi_ss.d(:,pos_p);
Psi.Dpi_w=Psi_ss.d(:,pos_w);
Psi.ss_pos=M_state_pos;
Psi.p_pos=ny1;
Psi.w_pos=ny2;


%%
function data_sgn_vec=sgn_vec(M_id,A_id,log,data_sgn_vec)

if log(M_id,4)~=27
    data_sgn_vec=[data_sgn_vec M_id];
    return
end

% data_sgn_vec=[data_sgn_vec M_id];

id1=log(M_id,6);
id2=log(M_id,7);

N1=0;
if id1==0
    N1=1;
else
    if log(id1,1)==1
        N1=1;
    else
        if log(id1,4)~=27
            N1=1;
            data_sgn_vec=[data_sgn_vec id1];
        end
    end
end

N2=0;
if id2==0
    N2=1;
else
    if log(id2,1)==1
        N2=1;
    else
        if log(id2,4)~=27
            N2=1;
            data_sgn_vec=[data_sgn_vec id2];
        end
    end
end

N3=0;
if A_id==0
    N3=1;
else
    if log(A_id,1)==1
        N3=1;
    else
        if log(A_id,4)~=27
            N3=1;
        end
    end
end

cond = mat2str([N1 N2 N3]);
switch cond
    case '[1 1 1]'
        return
    case '[1 1 0]'
        id1 = A_id;
        id2 = 0;
        A_id = 0;
    case '[1 0 1]'
        id1 = id2;
        id2 = 0;
        A_id = 0;
    case '[1 0 0]'
        id1 = id2;
        id2 = A_id;
        A_id = 0;
    case '[0 1 1]'
        id2 = 0;
        A_id = 0;
    case '[0 1 0]'
        id2 = A_id;
        A_id = 0;
    case '[0 0 1]'
        A_id = 0;
        %     case '[0 0 0]'
end
data_sgn_vec=sgn_vec(id1,id2,log,data_sgn_vec);