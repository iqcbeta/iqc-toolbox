%% now we define the "non-KYP" lmi
if vrb,
    disp_str(53)
end
str{sc}='\n%% Define non-KYP lmis ...';
sc=sc+1;

kvar=find(E.T==lmi);

switch ABST.lmitool
    case 'lmilab'
        for k=1:length(kvar),
            ks=num2str(k);
            kk=kvar(k);
            [nn,mm]=size(E.AB{kk});    % dimensions
            nns=num2str(nn);
            vs=mm-nn;                       % LMI size
            vss=num2str(vs);
            if nn>0,
                eval(['state' ks '=E.AB{kk};']);
                str{sc}=['p',ks,'=lmivar(1,[',nns,' 1]);']; %#ok<*AGROW>
                eval(str{sc});
                sc=sc+1;
                if strcmp(ABST.systemtype,'continuous'),
                    str{sc}=['lmiterm([',ks,' 1 1 p',ks,...
                        '],[eye(',nns,');zeros(',vss,',',nns,...
                        ')],state',ks,',''s'');'];
                    eval(str{sc});
                    sc=sc+1;
                elseif strcmp(ABST.systemtype,'discrete'),
                    str{sc}=['lmiterm([',ks,' 1 1 p',ks,...
                        '],state',ks,''',state',ks,');'];
                    eval(str{sc});
                    sc=sc+1;
                    str{sc}=['lmiterm([',ks,' 1 1 p',ks,...
                        '],[-eye(',nns,') zeros(',nns,',',vss,...
                        ')]'',[eye(',nns,') zeros(',nns,',',vss,...
                        ')]);'];
                    eval(str{sc});
                    sc=sc+1;
                else
                    disp_str(55)
                end
            end
            nh=0;
            ng=sum(E.X{kk}(3,:));
            for i=1:size(E.X{kk},2),
                is=num2str(i);
                nvs=num2str(E.X{kk}(1,i));   % variable number
                eval(['l_' ks '_' is '=E.C{kk}(ng+1:ng+E.X{kk}(2,i),:)'';'])
                eval(['r_' ks '_' is '=E.C{kk}(nh+1:nh+E.X{kk}(3,i),:);'])
                str{sc}=['lmiterm([',ks,' 1 1 ',nvs,'],l_',ks,'_',is,...
                    ',r_',ks,'_',is,',''s'');']; %#ok<AGROW>
                eval(str{sc});
                sc=sc+1;
                nh=nh+E.X{kk}(3,i);
                ng=ng+E.X{kk}(2,i);
            end
        end
    case 'yalmip'
        for k=1:length(kvar);
            ks=num2str(k);
            kk=kvar(k);
            [nn,mm]=size(E.AB{kk});    % dimensions
            nns=num2str(nn);
            mms=num2str(mm);
            vs=mm-nn;
            vss=num2str(vs);
            str{sc}=['LMI',ks,'=zeros(',mms,');'];
            eval(str{sc});
            sc=sc+1;
            if nn>0;
                eval(['state' ks '=E.AB{kk};']);
                str{sc}=['p',ks,'=sdpvar(',nns,',',nns,',''symmetric'');'];
                eval(str{sc});
                sc=sc+1;
                num_var_pos=[num_var_pos;{['p',ks],'1',[]}];
                str{sc}=['ouf',ks,'=[state',ks,...
                    ';eye(',nns,') zeros(',nns,',',vss,')];'];
                eval(str{sc});
                sc=sc+1;
                str{sc}=['LMI',ks,'=LMI',ks,...
                    '+ouf',ks,'''*kron(',psi,',p',ks,')*ouf',ks,';'];
                eval(str{sc});
                sc=sc+1;
            end
            nh=0;
            ng=sum(E.X{kk}(3,:));
            for i=1:size(E.X{kk},2),
                is=num2str(i);
                vr=E.X{kk}(1,i);
                eval(['l_' ks '_' is '=E.C{kk}(ng+1:ng+E.X{kk}(2,i),:)'';'])
                eval(['r_' ks '_' is '=E.C{kk}(nh+1:nh+E.X{kk}(3,i),:);'])
                if vr < 0,
                    vrs=num2str(-vr);
                    str{sc}=['LMI',ks,'=LMI',ks,'+l_',ks,'_',is,'*X',vrs,...
                        '''*r_',ks,'_',is,'+(l_',ks,'_',is,'*X',vrs,...
                        '''*r_',ks,'_',is,')'';'];
                    eval(str{sc});
                    sc=sc+1;
                else
                    vrs=num2str(vr);
                    str{sc}=['LMI',ks,'=LMI',ks,'+l_',ks,'_',is,'*X',vrs,...
                        '*r_',ks,'_',is,'+(l_',ks,'_',is,'*X',vrs,...
                        '*r_',ks,'_',is,')'';'];
                    eval(str{sc});
                    sc=sc+1;
                end
                nh=nh+E.X{kk}(3,i);
                ng=ng+E.X{kk}(2,i);
            end
            str{sc}=['ALL_LMI = ALL_LMI + [LMI',ks,'<= -',sft_v,'*eye(',...
                mms,')];'];
            eval(str{sc});
            sc=sc+1;
        end
end