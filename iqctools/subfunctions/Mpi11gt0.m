str{sc}='\n%% (#)*Mpi*Psi_q>0 terms...';
sc=sc+1;

switch ABST.lmitool
    case 'lmilab'
        
        if size(A1,1)~=0
            nns=num2str(size(A1,1));
            str{sc}=['[Xpsi,nXpsi,sXpsi]=lmivar(1,[',nns,' 1]);'];
            eval(str{sc});
            sc=sc+1;
            
            vss=num2str(size(B1,2));
            
            
            switch ABST.systemtype
                case 'continuous'
                    str{sc}=['lmiterm([-',ks,' 1 1 Xpsi],1,A1,''s'');'];
                    eval(str{sc});
                    sc=sc+1;
                    
                    str{sc}=['lmiterm([-',ks,' 1 2 Xpsi],1,B1);'];
                    eval(str{sc});
                    sc=sc+1;
                    
                case 'discrete'
                    str{sc}=['lmiterm([-',ks,' 1 1 Xpsi],A1'',A1);'];
                    eval(str{sc});
                    sc=sc+1;
                    
                    str{sc}=['lmiterm([-',ks,' 1 1 Xpsi],-1,1);'];
                    eval(str{sc});
                    sc=sc+1;
                    
                    str{sc}=['lmiterm([-',ks,' 1 2 Xpsi],A1'',B1);'];
                    eval(str{sc});
                    sc=sc+1;
                    
                    str{sc}=['lmiterm([-',ks,' 2 2 Xpsi],B1'',B1);'];
                    eval(str{sc});
                    sc=sc+1;
            end
        end
        
%         %%
%         str{sc}=['LMI',ks,'_11=zeros(',num2str(size(A1,1)),');'];
%         sc=sc+1;
%         str{sc}=['LMI',ks,'_12=zeros(',num2str(size(A1,1)),',',num2str(size(B1,2)),');'];
%         sc=sc+1;
%         str{sc}=['LMI',ks,'_22=zeros(',num2str(size(B1,2)),');'];
%         sc=sc+1;
        %%
        
        % Mpi
        for j=1:E.nsimple,
            js=num2str(j);
            nh=0;
            ng=0;
            n=E.simples(j);
            for i=1:size(E.X{n},2),
                is=num2str(i);
                vrs = num2str(E.X{n}(1,i));
                if size(A1,1)~=0
                    str{sc}=['lmiterm([-',ks,' 1 1 ',vrs,...
                        '],C1''*MpiL_',js,'_',is,...
                        ',MpiR_',js,'_',is,'*C1,''s'');'];
                    eval(str{sc})
                    sc=sc+1;
                    
                    str{sc}=['lmiterm([-',ks,' 1 2 ',vrs,...
                        '],C1''*MpiL_',js,'_',is,...
                        ',MpiR_',js,'_',is,'*D1);'];
                    eval(str{sc})
                    sc=sc+1;
                    
                    str{sc}=['lmiterm([-',ks,' 1 2 -',vrs,...
                        '],C1''*MpiR_',js,'_',is,...
                        ''',MpiL_',js,'_',is,'''*D1);'];
                    eval(str{sc})
                    sc=sc+1;
                    
                    str{sc}=['lmiterm([-',ks,' 2 2 ',vrs,...
                        '],D1''*MpiL_',js,'_',is,...
                        ',MpiR_',js,'_',is,'*D1,''s'');'];
                    eval(str{sc})
                    sc=sc+1;
                    
                else
                    str{sc}=['lmiterm([-',ks,' 1 1 ',vrs,...
                        '],D1''*MpiL_',js,'_',is,...
                        ',MpiR_',js,'_',is,'*D1,''s'');'];
                    eval(str{sc})
                    sc=sc+1;
                end
                
                %%
                %                 if E.X{n}(1,i)>0
                %                     str{sc}=['LMI',ks,'_11=LMI',ks,'_11+(C1''*MpiL_',js,'_',is,...
                %                         '*X',vrs,'*MpiR_',js,'_',is,'*C1)+(C1''*MpiL_',js,'_',is,...
                %                         '*X',vrs,'*MpiR_',js,'_',is,'*C1)'''];
                %                     sc=sc+1;
                %
                %                     str{sc}=['LMI',ks,'_12=LMI',ks,'_12+(C1''*MpiL_',js,'_',is,...
                %                         '*X',vrs,'*MpiR_',js,'_',is,'*D1)+(C1''*MpiR_',js,'_',is,...
                %                         '*X',vrs,'''*MpiL_',js,'_',is,'''*D1)'];
                %                     sc=sc+1;
                %
                %                     str{sc}=['LMI',ks,'_22=LMI',ks,'_22+(D1''*MpiL_',js,'_',is,...
                %                         '*X',vrs,'*MpiR_',js,'_',is,'*D1)+(D1''*MpiL_',js,'_',is,...
                %                         '*X',vrs,'*MpiR_',js,'_',is,'*D1)'''];
                %                     sc=sc+1;
                %                 else
                %                     vrs = num2str(-E.X{n}(1,i));
                %                     str{sc}=['LMI',ks,'_11=LMI',ks,'_11+(C1''*MpiL_',js,'_',is,...
                %                         '*X',vrs,'''*MpiR_',js,'_',is,'*C1)+(C1''*MpiL_',js,'_',is,...
                %                         '*X',vrs,'''*MpiR_',js,'_',is,'*C1)'''];
                %                     sc=sc+1;
                %
                %                     str{sc}=['LMI',ks,'_12=LMI',ks,'_12+(C1''*MpiL_',js,'_',is,...
                %                         '*X',vrs,'''*MpiR_',js,'_',is,'*D1)+(C1''*MpiR_',js,'_',is,...
                %                         '''*X',vrs,'*MpiL_',js,'_',is,'''*D1)'];
                %                     sc=sc+1;
                %
                %                     str{sc}=['LMI',ks,'_22=LMI',ks,'_22+(D1''*MpiL_',js,'_',is,...
                %                         '*X',vrs,'''*MpiR_',js,'_',is,'*D1)+(D1''*MpiL_',js,'_',is,...
                %                         '*X',vrs,'''*MpiR_',js,'_',is,'*D1)'''];
                %                     sc=sc+1;
                %                 end
                %%
                
                
                nh=nh+E.X{n}(3,i);
                ng=ng+E.X{n}(2,i);
            end
        end
    case 'yalmip'
        nns=num2str(size(A1,1));
        str{sc}=['Xpsi=sdpvar(',nns,');'];
        eval(str{sc});
        sc=sc+1;
        num_var_pos=[num_var_pos;{'Xpsi','1',[]}];
        
        vss=num2str(size(B1,2));
        str{sc}=['ouf_Psi_q=[eye(',nns,') zeros(',nns,',',vss,');'...
            'A1 B1;C1 D1];'];
        eval(str{sc});
        sc=sc+1;
        
        mms=num2str(size(Mpi,1));
        switch ABST.systemtype
            case 'continuous'
                str{sc}=['M_Xpsi_Mpi=[zeros(',nns,',',nns,...
                    ') Xpsi zeros(',nns,',',mms,');'...
                    'Xpsi zeros(',nns,',',nns,') zeros(',nns,',',mms,');'...
                    'zeros(',mms,',',cal_str(nns,nns),') Mpi];'];
            case 'discrete'
                str{sc}=['M_Xpsi_Mpi=[-Xpsi zeros(',nns,...
                    ',',cal_str(nns,mms),');'...
                    'zeros(',nns,',',nns,') Xpsi zeros(',nns,',',mms,');'...
                    'zeros(',mms,',',cal_str(nns,nns),') Mpi];'];
        end
        eval(str{sc});
        sc=sc+1;
        
        str{sc}='Mpi11=ouf_Psi_q''*M_Xpsi_Mpi*ouf_Psi_q;';
        eval(str{sc});
        sc=sc+1;
        
        if isa(Mpi11,'sdpvar') 
            str{sc}=['ALL_LMI=ALL_LMI+[Mpi11>=',sft_v,'*eye(',...
                num2str(size(Mpi11,1)),')];'];
            eval(str{sc});
            sc=sc+1;
        end
end
