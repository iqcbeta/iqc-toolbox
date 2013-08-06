if vrb,
    disp_str(52)
end
str{6}='\n%% Define multiplier variables ...';
sc=7; % string counter

kvar=find(E.T==var);

switch ABST.lmitool
    case 'lmilab'
        for k=1:length(kvar);
            ks=num2str(k);
            eval(['struct' ks '=E.L{' num2str(kvar(k)) '};'])
            str{sc}=['[X',ks,',nX',ks,',sX',ks,...
                ']=lmivar(3,struct' ks ');'];
            eval(str{sc});
            sc=sc+1;
        end
    case 'yalmip'
        vc=0;
        for k=1:length(kvar),
            ks   = num2str(k);
            varn = num2str(E.X{kvar(k)}(2,:));
            varm = num2str(E.X{kvar(k)}(3,:));
            switch ABST.log(kvar(k),5)
                case 1  % symmetric
                    str{sc}=['X' ks '=sdpvar(' varn ',' varm ...
                        ',''symmetric'');']; %#ok<*AGROW>
                    eval(str{sc});
                    sc=sc+1;
                    num_var_pos=[num_var_pos;{['X',ks],'1',[]}];
                case 2  % rectangular
                    str{sc}=['X' ks '=sdpvar(' varn ',' varm ',''full'');'];
                    eval(str{sc});
                    sc=sc+1;
                    num_var_pos=[num_var_pos;{['X',ks],'2',[]}];
                case 3  % diagonal
                    str{sc} = ['X' ks '=diag(sdpvar(' varn ',1));'];
                    eval(str{sc});
                    sc=sc+1;
                    num_var_pos=[num_var_pos;{['X',ks],'3',[]}];
                case 4  % skew symmetric
                    str{sc}=['X' ks '=sdpvar(' varn ',' varm ',''skew'');'];
                    eval(str{sc});
                    sc=sc+1;
                    num_var_pos=[num_var_pos;{['X',ks],'4',[]}];
                case 5  % variable
                    str{sc} = ['X' ks '=zeros(' varn ',' varm ');'];
                    eval(str{sc});
                    sc=sc+1;
                    eval(['struct' ks '= E.L{' num2str(kvar(k)) '};']);
                    M   = E.L{kvar(k)};
                    Ma  = sort(M(:));
                    Ms  = Ma(Ma>0);
                    Msc = Ms(1);
                    cnt = 1;
                    for flcnt = 2: length(Ms)
                        if Ms(flcnt) > Msc(cnt);
                            Msc(cnt+1) = Ms(flcnt);
                            cnt = cnt + 1;
                        end
                    end
                    allvc=[];
                    for flcnt = 1: length(Msc)
                        Mscvar = num2str(Msc(flcnt));
                        vc=vc+1;
                        vcs = num2str(vc);
                        str{sc}  = ['y' vcs '= sdpvar;'];
                        eval(str{sc});
                        sc=sc+1;
                        str{sc} = ['X' ks '= X' ks '+((abs(struct' ks...
                            ')==' Mscvar ').*sign(struct' ks '))*y' vcs...
                            ';'];
                        eval(str{sc});
                        sc=sc+1;
                        allvc=[allvc;['y',vcs]];
                    end
                    num_var_pos=[num_var_pos;{['X',ks],'5',...
                        num2str(allvc)}];
                otherwise
                    switch ABST.log(kvar(k),4),
                        case 28, % subsref
                            rc=find(kvar==ABST.log(kvar(k),6));
                            p_rc=find(E.L{kvar(rc)}==E.L{kvar(k)});
                            str{sc} = ['X' ks '= X' num2str(rc)...
                                '(' num2str(p_rc(1)) ');' ];
                            eval(str{sc});
                            sc=sc+1;
                        case 29, % subsasgn
                            rc=find(kvar==ABST.log(kvar(k),7));
                            p_rc=find(E.L{kvar(rc)}==E.L{kvar(k)});
                            str{sc} = ['X' ks '= X' num2str(rc)...
                                '(' num2str(p_rc(1)) ');' ];
                            eval(str{sc});
                            sc=sc+1;
                        otherwise
                            disp_str(55)
                    end
            end
        end   
end