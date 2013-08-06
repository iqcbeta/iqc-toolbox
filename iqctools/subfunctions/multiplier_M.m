str{sc}='\n%% define Mpi terms...';
sc=sc+1;

%% Multipliers: Mpi
% [ 0  I   0    0 ]^T      [ 0  I   0    0 ]
% [ Cq 0  Dqw  Dqp]   Mpi  [ Cq 0  Dqw  Dqp] = L_j_i*X*R_j_i = IQCpi
% [ 0  0   I    0 ]        [ 0  0   I    0 ]
% [ 0  0   0    I ]        [ 0  0   0    I ]

nns=num2str(size(Cpi,1));

str{sc}=['Mpi=zeros(',nns,',',nns,');'];
eval(str{sc});
sc=sc+1;

%% main_ouf
%            [ 0  I   0    0 ]
% main_ouf = [ Cq 0  Dqw  Dqp]
%            [ 0  0   I    0 ]
%            [ 0  0   0    I ]

ny1=size(Api,1);
ny2=size(Cq,1);
ny3=size(Bw,2);
ny4=size(Bp,2);
ny=ny1+ny2+ny3+ny4;

main_ouf=zeros(ny,size(E.ab,2));
main_ouf(1:ny1,Pi_pos)=eye(ny1);
main_ouf(ny1+1:ny1+ny2,N_pos)=Cq;
main_ouf(ny1+1:ny1+ny2,w_pos)=Dqw;
main_ouf(ny1+1:ny1+ny2,p_pos)=Dqp;
main_ouf(ny1+ny2+1:ny1+ny2+ny3,w_pos)=eye(ny3);
main_ouf(ny1+ny2+ny3+1:end,p_pos)=eye(ny4);


%% MpiR_*_* MpiL_*_*
for j=1:E.nsimple,
    js=num2str(j);
    nh=0;
    ng=0;
    n=E.simples(j);
    for i=1:size(E.X{n},2),
        is=num2str(i);
        
        eval(['L_' js '_' is ...
            '=E.B{E.simples(j)}(ng+1:ng+E.X{E.simples(j)}(2,i),:)'';'])
        eval(['R_' js '_' is ...
            '=E.gqfm(j)*E.C{E.simples(j)}(nh+1:nh+' ...
            'E.X{E.simples(j)}(3,i),:);'])
        
        eval(['MpiL_',js,'_',is,'=(main_ouf''\L_',js,'_',is,');']);
        eval(['MpiR_',js,'_',is,'=(R_',js,'_',is,'/main_ouf);']);
        
        if ~isempty(find(eval(['L_',js,'_',is,])-...
                main_ouf'*eval(['MpiL_',js,'_',is])>eps, 1))
            disp_str(72)
        end
        
        if ~isempty(find(eval(['R_',js,'_',is,])-...
                eval(['MpiR_',js,'_',is])*main_ouf>eps, 1))
            disp_str(72)
        end
        
        nh=nh+E.X{n}(3,i);
        ng=ng+E.X{n}(2,i);
    end
end
