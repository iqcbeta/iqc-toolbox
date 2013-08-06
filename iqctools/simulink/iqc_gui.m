function [gain,exe_comment]=iqc_gui(sys_name,datafile_name,datafile_type)

old_path=pwd;

global ext_iqc_setvalues
if isempty(ext_iqc_setvalues)
    
    if isempty(strfind(sys_name,'.mdl'))
        OriginallyFileName_full=[sys_name,'.mdl'];
        
    else
        OriginallyFileName_full=sys_name
    end
    file_path=which(OriginallyFileName_full,'-all');
    if size(file_path,1)~=1
        choose_n=0;
        choose_all=size(file_path,1);
        disp_str(59,'the same name in different directories are found')
        for i1=1:choose_all
            disp([sft(num2str(i1),4),':',file_path{i1}])
        end
        while choose_n==0
            choose_n=input(['(1-',num2str(choose_all),') ?: '], 's');
            if isempty(find(1:choose_all==str2double(choose_n), 1))
                choose_n=0;
            else
                choose_n=str2double(choose_n);
            end
        end
        file_path=file_path{choose_n};
    else
        file_path=file_path{1};
    end
    dot_pos=(strfind(file_path,OriginallyFileName_full));
    cd(file_path(1:dot_pos-1));
    
    sim(sys_name)
    return
end

if nargin==0,
    disp('You have to build a SIMULINK block diagram')
    disp('using the blocks from "iqc_lib" before')
    disp('using this program. ')
    sys_name=input('Enter diagram "filename" now: ','s');
    disp(' ')
elseif nargin==1,
    dot_pos=strfind(sys_name,'.mdl');
    if ~isempty(dot_pos),
        sys_name=sys_name(1:dot_pos(1,1)-1);
    end
    datafile_name=[];
    isdatafile=0;
elseif nargin==2,
    isdatafile=1;
    dot_pos=strfind(datafile_name,'.m');
    if ~isempty(dot_pos)
        datafile_name=datafile_name(1:dot_pos(1,1)-1);
    end
    datafile_name_full=[datafile_name,'.m'];
elseif nargin==3,
    isdatafile=1;
    dot_pos=(strfind(datafile_name,'.'));
    if ~isempty(dot_pos)
        datafile_name=datafile_name(1:dot_pos(1,1)-1);
    end
    switch datafile_type,
        case 1,
            datafile_name_full=[datafile_name,'.m'];
        case 2,
            datafile_name_full=[datafile_name,'.mat'];
        otherwise
            error('Data file can be only MATLAB M-file or Mat-datafile')
    end
else
    error('Too many inputs')
end

if isdatafile,
    if ~exist(datafile_name_full,'file')
        disp('The specified data file does not exist ...')
        error('Or the file is not in the directory which is on the MATLAB search path')
    end
end

if isdatafile,
    keep_linfile=1;
else
    keep_linfile=0;
end

[exe_comment,ecctr,LinearFileName,save_lin_com]=...
    gui_iqctbx(sys_name,keep_linfile);

exe_sys_name=[sys_name,'_iqcexe'];
exe_sys_name_full=[sys_name,'_iqcexe.m'];

fp=fopen(exe_sys_name_full,'wt');
fprintf(fp,['function gain=',exe_sys_name,'\n']);

if ~isdatafile
    for forlp_counter=1:ecctr
        fprintf(fp,[exe_comment{forlp_counter},'\n']);
    end
else
    if type==1
        fp1=fopen(datafile_name_full,'rt');
        while ~feof(fp1),
            line_temp=fgets(fp1);
            fprintf(fp,[line_temp,'\n']);
        end
        fclose(fp1);
        for forlp_counter=1:ecctr
            fprintf(fp,[exe_comment{forlp_counter},'\n']);
        end
    elseif type==2
        fprintf(fp,['load ',datafile_name,'\n']);
        for forlp_counter=1:ecctr
            fprintf(fp,[exe_comment{forlp_counter},'\n']);
        end
    end
end
fclose(fp);

if isdatafile==0,
    
    for forlp_counter=1:ecctr
        eval(exe_comment{forlp_counter})
    end
    if ~save_lin_com
        delete(exe_sys_name_full);
        delete([LinearFileName,'.mdl']);
    end
else
    disp(' ')
    strtmp1='iqc_gui.m has processed the graphic file ';
    strtmp2='... and produced the corresponding executable M-file ';
    disp([strtmp1,sys_name,'.mdl'])
    disp([strtmp2,exe_sys_name_full])
end

old_sim_name={'r14','r14sp1','r14sp2','r14sp3','r2006a','r2006b',...
    'r2007a','r2007b','r2008a','r2008b','r2009a','r2009b',...
    'r2010a','r2010b','r2011a','r2011b','r2012a','r2012b','unknown_version'};
for i1=1:length(old_sim_name)
    if exist([sys_name,'_lin.mdl.',old_sim_name{i1}])
        delete([sys_name,'_lin.mdl.',old_sim_name{i1}]);
    end
    if exist([sys_name,'.mdl.',old_sim_name{i1}])
        delete([sys_name,'.mdl.',old_sim_name{i1}]);
    end
end
cd(old_path);

%%
function [exe_comment,sc,LinearFileName,save_lin_com]=gui_iqctbx(filename,keep_linfile) %#ok<INUSD>
%%
if nargin==1,
    keep_linfile=0;
end

exe_comment={};        % executable comments cell array
sc=0;                  % executable comments cell array counter

disp(' ')
disp(' ')
disp('*************************************************')
disp(['iqc_gui: processing diagram "' filename '.mdl" ...'])
disp('*************************************************')
disp(' ')

OriginallyFileName_full=[filename '.mdl'];
LinearFileName=[filename '_lin'];
LinearFileName_full=[LinearFileName '.mdl'];
if exist(LinearFileName_full,'file'),
    delete(LinearFileName_full);
end

%-------------------------------copy file and rename file_lin
file_path=which(OriginallyFileName_full,'-all');
file_path=file_path{1};


dot_pos=(strfind(file_path,OriginallyFileName_full));
cd(file_path(1:dot_pos-1));
copyfile(OriginallyFileName_full,LinearFileName_full)

load_system(LinearFileName_full);

IQC_INFO=[];      % containing information of IQC block
iqc_block=[];
A=zeros(1,12);

global ext_iqc_setvalues
iqc_block_name=fieldnames(ext_iqc_setvalues);

% 用於 多個 performance 方塊使用，暫無效果
for i1=1:length(iqc_block_name),
    if strcmp(iqc_block_name{i1},'runguiiqc')
        iqc_block_name(i1)='';
        break;
    end
end
%---------------------------------------

name=[];
position=[];

counter=0;
sampletimeset=[];

exp_block_name={'ext_sim_saturation','ext_sim_deadzone',...
    'ext_sim_intenc_deadzone','ext_sim_uncertain_delay','ext_sim_normboundsystem',...
    'ext_sim_unknown_gain','ext_sim_rate_limiter','ext_sim_uncertain_cdelay'};
white_number=0;

L2_block=find_system(LinearFileName,'LookUnderMasks','all',...
    'FunctionName','ext_iqc_L2signal');
White_block=find_system(LinearFileName,'LookUnderMasks','all',...
    'FunctionName','ext_iqc_white');

new_block_item='';
new_block=[];
for i1=1:length(White_block)
    new_block_item=White_block{i1,1};
    dot_pos=strfind(new_block_item,'/');
    new_block_item=new_block_item(1:dot_pos(length(dot_pos))-1);
    new_block=[new_block;{new_block_item}];
end
White_block=new_block;
clear new_block_item dot_pos

Signal_block=[L2_block;White_block];
L2block_name=get_param(Signal_block,'name');

white_pos=[];
L2signal_pos=[];
performance_pos=[];


for i1=1:length(iqc_block_name)
    idx1=1;
    iqc_parameter_name='';
    new_iqc_block=[];
    
    if strcmp(iqc_block_name{i1},'ext_iqc_white')
        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.idx'];
        for i2=1:eval(str)
            str=['ext_iqc_setvalues.',iqc_block_name{i1},'.newpath{i2}'];
            new_iqc_block=eval(str);
            iqc_block=[iqc_block;{new_iqc_block}];
            str=['ext_iqc_setvalues.',iqc_block_name{i1},'.sampletime{i2}(2)'];
            sampletimeset=[sampletimeset eval(str)];
            str=['ext_iqc_setvalues.',iqc_block_name{i1}];
            iqc_parameter_name=fieldnames(eval(str));
            str=['ext_iqc_setvalues.',iqc_block_name{i1},'.name{i2}'];
            IQC.name{counter+i2,1}=eval(str);
            IQC.type{counter+i2,1}=iqc_block_name{i1};
            IQC.rela{counter+i2,2}=[0 0];
            white_number=[white_number i1];
            
            str=['ext_iqc_setvalues.',iqc_block_name{i1},'.newpath{i2}'];
            A(1,1:4)=get_param(eval(str),'Position');
            str=['ext_iqc_setvalues.',iqc_block_name{i1},'.sig_size{1}(2)'];
            A(1,6)=eval(str);
            A(1,5)=0;
            IQC_INFO=[IQC_INFO;A];
            white_pos=[white_pos size(IQC_INFO,1)];
            idx2=1;
            parameter='';
            for i3=1:length(iqc_parameter_name), %path,name,size,sampletime,idx
                switch iqc_parameter_name{i3}
                    case 'path'
                    case 'name'
                    case 'sampletime'
                    case 'idx'
                    otherwise
                        parameter.name{idx2}=iqc_parameter_name{i3};
                        str=['ext_iqc_setvalues.',iqc_block_name{i1},...
                            '.',iqc_parameter_name{i3},'(i2)'];
                        parameter.value{idx2}=eval(str);
                        idx2=idx2+1;
                end
            end
            info=IQC_par(counter+i2,iqc_block_name{i1},parameter);
            IQC.rela{counter+i2,1}=info(1,:);
        end
        counter=counter+i2;
        continue
    end
    
    if strcmp(iqc_block_name{i1},'ext_iqc_performance')
        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.newpath{1}'];
        new_iqc_block=eval(str);
        iqc_block=[iqc_block;{new_iqc_block}];
        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.idx'];
        if eval(str)>1
            error('performance block error!')
        end
        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.sampletime{1}(1)'];
        %         sampletimeset=[sampletimeset eval(str)];
        str=['ext_iqc_setvalues.',iqc_block_name{i1}];
        iqc_parameter_name=fieldnames(eval(str));
        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.newpath{1}'];
        A(1,1:4)=get_param(eval(str),'Position');
        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.size{1}(1)'];
        A(1,5)=eval(str);
        A(1,6)=0;
        IQC_INFO=[IQC_INFO;A];
        performance_pos=[performance_pos size(IQC_INFO,1)];
        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.name{1}'];
        IQC.name{counter+1,1}=eval(str);
        IQC.type{counter+1,1}=iqc_block_name{i1};
        IQC.rela{counter+1,2}=A(1,[6 5]);
        idx2=1;
        parameter='';
        for i3=1:length(iqc_parameter_name), %path,name,size,sampletime,idx
            switch iqc_parameter_name{i3}
                case 'path'
                case 'newpath'
                case 'name'
                case 'size'
                case 'jump_idx'
                case 'sampletime'
                case 'save_lin_com'
                    switch eval(['ext_iqc_setvalues.',iqc_block_name{i1},...
                            '.',iqc_parameter_name{i3},'{1}'])
                        case 0
                            save_lin_com=0;
                        case 1
                            save_lin_com=1;
                    end
                case 'idx'
                otherwise
                    parameter.name{idx2}=iqc_parameter_name{i3};
                    parameter.value{idx2}=eval(['ext_iqc_setvalues.',iqc_block_name{i1},...
                        '.',iqc_parameter_name{i3}]);
                    idx2=idx2+1;
            end
        end
        info=IQC_par(counter+1,iqc_block_name{i1},parameter);
        if size(info,1)~=1,
            IQC.lmi{counter+1,1}=info{2,:};
            IQC.rela{counter+1,1}=info{1,:};
        else
            IQC.rela{counter+1,1}=info(1,:);
        end
        counter=counter+1;
        continue
    end
    
    % others
    str=['ext_iqc_setvalues.',iqc_block_name{i1},'.idx'];
    for i2=1:eval(str)
        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.newpath{i2}'];
        new_iqc_block{idx1,1}=eval(str);
        idx1=idx1+1;
        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.sampletime{i2}(1:2)'];
        if strcmp(iqc_block_name{i1},'ext_iqc_L2signal')
            sampletimeset=[sampletimeset eval(str)];
        end
        str=['ext_iqc_setvalues.',iqc_block_name{i1}];
        iqc_parameter_name=fieldnames(eval(str));
        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.newpath{i2}'];
        A(1,1:4)=get_param(eval(str),'Position');
        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.size{i2}(1)'];
        A(1,5)=eval(str);
        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.size{i2}(2)'];
        A(1,6)=eval(str);
        IQC_INFO=[IQC_INFO;A];
        if strcmp(iqc_block_name{i1},'ext_iqc_L2signal'),
            L2signal_pos=[L2signal_pos size(IQC_INFO,1)];
        end
        parameter='';
        if strcmp(iqc_block_name{i1},'ext_iqc_custom'),
            % 計算方塊之自訂方塊
            str='ext_iqc_setvalues.ext_iqc_custom.list_iqc_rel{i2}';
            idx2=strfind(eval(str),'|');
            if isempty(idx2)
                IQC.name{counter+i2,1}=[char(get_param(new_iqc_block(i2),...
                    'Name')),' name1'];
                IQC.type{counter+i2,1}=iqc_block_name{i1};
                IQC.rela{counter+i2,2}=A(1,[6 5]);
                str='ext_iqc_setvalues.ext_iqc_custom.list_iqc_rel{i2}';
                IQC.rela{counter+i2,1}=eval(str);
            else
                idx3=1;
                IQC.name{counter+i2,1}=char(get_param(new_iqc_block(i2),...
                    'Name'));
                IQC.type{counter+i2,1}=iqc_block_name{i1};
                for i3=1:length(idx2)+1
                    IQC.rela{counter+i2,2}{i3}=A(1,[6 5]);
                    if i3~=length(idx2)+1
                        str=['ext_iqc_setvalues.ext_iqc_custom.',...
                            'list_iqc_rel{i2}(idx3:idx2(i3)-1)'];
                        IQC.rela{counter+i2,1}{i3}=eval(str);
                        idx3=idx2(i3)+1;
                    else
                        str=['ext_iqc_setvalues.ext_iqc_custom.',...
                            'list_iqc_rel{i2}'];
                        str=[str,'(idx3:size(',str,',2))'];
                        IQC.rela{counter+i2,1}{i3}=eval(str);
                    end
                end
            end
        elseif ~isempty(findcell(exp_block_name,iqc_block_name{i1}))
            % 模擬方塊
            str=['ext_iqc_setvalues.',iqc_block_name{i1},...
                '.iqctype{i2}'];
            idx2=length(eval(str));
            IQC.name{counter+i2,1}=char(get_param(new_iqc_block(i2),...
                'Name'));
            IQC.type{counter+i2,1}=iqc_block_name{i1};
            counter2=0;
            for i3=1:idx2
                str=['ext_iqc_setvalues.',iqc_block_name{i1},...
                    '.iqctype{i2}{i3}{1,1}'];
                type_name=eval(str);
                dot_pos=strfind(eval(str),'.');
                type_name=['ext_',type_name(1:dot_pos(1)-1)];
                if strcmp(type_name,'ext_iqc_custom')
                    str=['ext_iqc_setvalues.',iqc_block_name{i1},...
                        '.iqctype{i2}{i3}{1,2}'];
                    list_iqc_rel=eval(str);
                    idx3=strfind(list_iqc_rel,'|');
                    if isempty(idx3)
                        IQC.rela{counter+i2,2}{counter2+i3}=A(1,[6 5]);
                        IQC.rela{counter+i2,1}{counter2+i3}=list_iqc_rel;
                    else
                        idx4=1;
                        for i4=1:length(idx3)+1
                            IQC.rela{counter+i2,2}{counter2+i3+i4-1}=...
                                A(1,[6 5]);
                            if i4~=length(idx3)+1
                                IQC.rela{counter+i2,1}{counter2+i3+i4-1}=...
                                    list_iqc_rel(idx4:idx3(i4)-1);
                                idx4=idx3(i4)+1;
                            else
                                IQC.rela{counter+i2,1}{counter2+i3+i4-1}=...
                                    list_iqc_rel(idx4:end);
                            end
                        end
                        counter2=counter2+length(idx2);
                    end
                else
                    str=['ext_iqc_setvalues.',iqc_block_name{i1},...
                        '.iqctype{i2}{i3}'];
                    idx2=1;
                    for i4=1:size(eval(str),1)
                        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.iqctype{i2}{i3}{i4,1}'];
                        dot_pos=strfind(eval(str),'.');
                        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.iqctype{i2}{i3}{i4,1}'];
                        str=[str,'(dot_pos+1:size(',str,',2))'];
                        parameter.name{idx2}=eval(str);
                        str=['ext_iqc_setvalues.',iqc_block_name{i1},'.iqctype{i2}{i3}(i4,2)'];
                        parameter.value{idx2}=eval(str);
                        idx2=idx2+1;
                    end
                    info=IQC_par(counter+i2,type_name,parameter);
                    IQC.rela{counter+i2,2}{counter2+i3}=A(1,[6 5]);
                    IQC.rela{counter+i2,1}{counter2+i3}=info(1,:);
                end
            end
        else %計算方塊之其它
            IQC.name{counter+i2,1}=char(get_param(new_iqc_block(i2),'Name'));
            IQC.type{counter+i2,1}=iqc_block_name{i1};
            IQC.rela{counter+i2,2}=A(1,[6 5]);
            idx2=1;
            for i3=1:length(iqc_parameter_name), %path,name,size,sampletime,idx
                switch iqc_parameter_name{i3}
                    case 'path'
                    case 'newpath'
                    case 'name'
                    case 'size'
                    case 'sampletime'
                    case 'idx';
                    otherwise
                        parameter.name{idx2}=iqc_parameter_name{i3};
                        parameter.value{idx2}=eval(['ext_iqc_setvalues.',...
                            iqc_block_name{i1},'.',...
                            iqc_parameter_name{i3},'(i2)']);
                        idx2=idx2+1;
                end
            end
            if isempty(parameter)
                continue
            end
            info=IQC_par(counter+i2,iqc_block_name{i1},parameter);
            if ~isempty(info)
                IQC.rela{counter+i2,1}=info(1,:);
            else
                IQC.rela{counter+i2,1}='';
            end
        end
    end
    iqc_block=[iqc_block;new_iqc_block];
    counter=counter+i2;
end

% IQC block can not be in the Subsystem
for i1=1:length(iqc_block)
    dot_pos=strfind(iqc_block{i1},'/');
    if length(dot_pos)>2
        error('IQC block can not be in the subsystem')
        iqctool('simclear')
    end
end

% delete main floor inport and outport
main_in=find_system(LinearFileName,'SearchDepth',1,'BlockType','Inport');
delete_block(main_in);
main_out=find_system(LinearFileName,'SearchDepth',1,'BlockType','Outport');
delete_block(main_out);
save_system(LinearFileName,[],'OverwriteIfChangedOnDisk',true);

% additive Inport and Outport
IQC_INFO=add_blk(IQC_INFO,LinearFileName);
save_system(LinearFileName,[],'OverwriteIfChangedOnDisk',true);

% change the connect line
linechang(LinearFileName,iqc_block,IQC,IQC_INFO);
save_system(LinearFileName,[],'OverwriteIfChangedOnDisk',true);
close_system(LinearFileName);

IQC_INFO=getportdim(IQC_INFO,LinearFileName);

if length(unique(sampletimeset))~=1
    error('the block sample time does not match')
end

% Second Round:
% Defind the signals and corresponding IQC relationships of
% the signals, and run the solver
% if (option(5)==0)|(option(5)==777),
%     disp('define the system and run solver ...')
% end
sc=1;

if unique(sampletimeset)==0
    exe_comment{sc}='abst_init_iqc;';
else
    exe_comment{sc}='abst_init_iqc(1);';
end
sc=sc+1;

exe_comment{sc}=['load_system(''',LinearFileName,''');'];
sc=sc+1;
if unique(sampletimeset)==0
    exe_comment{sc}=['[A_iqcexe,B_iqcexe,C_iqcexe,D_iqcexe]=linmod2(''',...
        LinearFileName,''');'];
else
    exe_comment{sc}=['[A_iqcexe,B_iqcexe,C_iqcexe,D_iqcexe]=dlinmod(''',...
        LinearFileName,''',',num2str(unique(sampletimeset)),');'];
end
sc=sc+1;
exe_comment{sc}=['close_system(''',LinearFileName,''',1);'];
sc=sc+1;

% set performance out signal pos and lmi solver set
if isfield(IQC,'lmi')
    if ~isempty(IQC.lmi{performance_pos,1})
        exe_comment{sc}=IQC.lmi{performance_pos,1};
        sc=sc+1;
    end
end
fc1=num2str(IQC_INFO(performance_pos,9));
fc2=num2str(IQC_INFO(performance_pos,10));
OUT=['w_iqcexe(',fc1,':',fc2,')'];

% set in signal pos

IN_pos=[];
for i1=1:length(white_pos)
    fc1=num2str(IQC_INFO(white_pos(i1),11));
    fc2=num2str(IQC_INFO(white_pos(i1),12));
    IN_pos=[IN_pos eval([fc1,':',fc2])];
end
for i1=1:length(L2signal_pos)
    fc1=num2str(IQC_INFO(L2signal_pos(i1),11));
    fc2=num2str(IQC_INFO(L2signal_pos(i1),12));
    IN_pos=[IN_pos eval([fc1,':',fc2])];
end
IN=['f_iqcexe(',mat2str(IN_pos),')'];

str='f_iqcexe=';

aloutvariable='';
if ~isempty(white_pos)
    str=[str,'['];
    str1='';
    num1=0;
    for i1=1:size(IQC_INFO,1)
        if ~isempty(find(white_pos==i1))
            str=[str,str1];
            str1='';
            num1=0;
            rela=IQC.rela{i1,1};
            dot_pos=strfind(rela,'=');
            if ~strcmp(rela(dot_pos-1),']')
                white_signal_name=rela(1:dot_pos-1);
            else
                dot_pos1=strfind(rela(1:dot_pos-1),',');
                white_signal_name=rela(2:dot_pos1(1)-1);
                for i2=1:length(dot_pos1)
                    str2='';
                    if dot_pos1(i2)==dot_pos1(length(dot_pos1))
                        str2=['global ',rela(dot_pos1(i2)+1:dot_pos-2)];
                        aloutvariable=[aloutvariable;...
                            {rela(dot_pos1(i2)+1:dot_pos-2)}];
                    else
                        str2=['global ',rela(dot_pos1(i2)+1:dot_pos(i2+1)-1)];
                        aloutvariable=[aloutvariable;...
                            {rela(dot_pos1(i2)+1:dot_pos(i2+1)-1)}];
                    end
                    exe_comment{sc}=str2;
                    sc=sc+1;
                end
            end
            exe_comment{sc}=[rela,';'];
            sc=sc+1;
            str=[str,white_signal_name,';'];
        else
            if IQC_INFO(i1,6)~=0
                num1=num1+IQC_INFO(i1,6);
                str1=['signal(',num2str(num1),')',';'];
            else
                continue
            end
        end
    end
    str=[str,str1];
    str=[str(1:size(str,2)-1),'];'];
    exe_comment{sc}=str;
    sc=sc+1;
else
    exe_comment{sc}='fsize=size(B_iqcexe,2);';
    sc=sc+1;
    str=[str,'signal(fsize);'];
    exe_comment{sc}=str;
    sc=sc+1;
end

% LTI system outputs, whose components are inputs of iqc-blocks
if unique(sampletimeset)==0
    exe_comment{sc}=...
        'w_iqcexe=ss(A_iqcexe,B_iqcexe,C_iqcexe,D_iqcexe)*f_iqcexe;';
else
    exe_comment{sc}=['w_iqcexe=ss(A_iqcexe,B_iqcexe,C_iqcexe,D_iqcexe,',...
        num2str(unique(sampletimeset)),')*f_iqcexe;'];
end
sc=sc+1;

% --- Define IQC relationship ---
for i1=1:size(IQC_INFO,1)
    if ~isempty(find(unique([white_pos,L2signal_pos,performance_pos])==i1))
        continue
    else
        fc1=num2str(IQC_INFO(i1,9));
        fc2=num2str(IQC_INFO(i1,10));
        iqc_block_in=['w_iqcexe(',fc1,':',fc2,')'];
        fc1=num2str(IQC_INFO(i1,11));
        fc2=num2str(IQC_INFO(i1,12));
        iqc_block_out=['f_iqcexe(',fc1,':',fc2,')'];
        
        if iscell(IQC.rela{i1,1}) % 同時滿足多個 IQC
            for i2=1:length(IQC.rela{i1,1})
                rela=IQC.rela{i1,1}{i2};
                dot_pos=strfind(rela,'=');
                lef_com=rela(1:dot_pos(1)-1);
                rig_com=rela(dot_pos(length(dot_pos))+1:length(rela));
                
                dot_pos=strfind(rig_com,'(');
                dot_pos1=strfind(rig_com,',');
                if ~isempty(dot_pos1)
                    rig_com=[rig_com(1:dot_pos),...
                        iqc_block_in,rig_com(dot_pos1(1):length(rig_com))];
                else
                    rig_com=[rig_com(1:dot_pos),iqc_block_in,')'];
                end
                
                if strcmp(lef_com(length(lef_com)),']')
                    dot_pos=strfind(lef_com,',');
                    if ~isempty(dot_pos)
                        signal_name=lef_com(2:dot_pos(1)-1);
                        for i3=1:length(dot_pos)
                            if i3==length(dot_pos)
                                exe_comment{sc}=['global ',...
                                    lef_com(dot_pos(i3)+1:length(lef_com)-1)];
                                sc=sc+1;
                                aloutvariable=[aloutvariable;...
                                    {lef_com(dot_pos(i3)+1:length(lef_com)-1)}];
                            else
                                exe_comment{sc}=['global ',...
                                    lef_com(dot_pos(i3)+1:dot_pos(i3+1)-1)];
                                sc=sc+1;
                                aloutvariable=[aloutvariable;...
                                    {lef_com(dot_pos(i3)+1:dot_pos(i3+1)-1)}];
                            end
                        end
                    else
                        signal_name=lef_com(2:length(lef_com)-1);
                    end
                    str=[lef_com,'=',rig_com,';'];
                    exe_comment{sc}=str;
                    sc=sc+1;
                    exe_comment{sc}=[iqc_block_out,'==',signal_name,';'];
                    sc=sc+1;
                else
                    signal_name=lef_com;
                    str=[iqc_block_out,'==',rig_com,';'];
                    exe_comment{sc}=str;
                    sc=sc+1;
                end
                
            end
            
        else %只滿足一個 IQC
            rela=IQC.rela{i1,1};
            dot_pos=strfind(rela,'=');
            lef_com=rela(1:dot_pos(1)-1);
            rig_com=rela(dot_pos(length(dot_pos))+1:length(rela));
            
            dot_pos=strfind(rig_com,'(');
            dot_pos1=strfind(rig_com,',');
            if ~isempty(dot_pos1)
                rig_com=[rig_com(1:dot_pos),...
                    iqc_block_in,rig_com(dot_pos1(1):length(rig_com))];
            else
                rig_com=[rig_com(1:dot_pos),iqc_block_in,')'];
            end
            
            if strcmp(lef_com(length(lef_com)),']')
                dot_pos=strfind(lef_com,',');
                if ~isempty(dot_pos)
                    signal_name=lef_com(2:dot_pos(1)-1);
                    for i3=1:length(dot_pos)
                        if i3==length(dot_pos)
                            exe_comment{sc}=['global ',...
                                lef_com(dot_pos(i3)+1:length(lef_com)-1)];
                            sc=sc+1;
                            aloutvariable=[aloutvariable;...
                                {lef_com(dot_pos(i3)+1:length(lef_com)-1)}];
                        else
                            exe_comment{sc}=['global ',...
                                lef_com(dot_pos(i3)+1:dot_pos(i3+1)-1)];
                            sc=sc+1;
                            aloutvariable=[aloutvariable;...
                                {lef_com(dot_pos(i3)+1:dot_pos(i3+1)-1)}];
                        end
                    end
                else
                    signal_name=lef_com(2:length(lef_com)-1);
                end
                str=[lef_com,'=',rig_com,';'];
                exe_comment{sc}=str;
                sc=sc+1;
                exe_comment{sc}=[iqc_block_out,'==',signal_name,';'];
                sc=sc+1;
            else
                signal_name=lef_com;
                str=[iqc_block_out,'==',rig_com,';'];
                exe_comment{sc}=str;
                sc=sc+1;
            end
        end
    end
end

exe_comment{sc}=['gain=iqc_gain_tbx(',IN,',',OUT,')'];
sc=sc+1;

if ~isempty(aloutvariable)
    for i1=1:length(aloutvariable)
        str=[aloutvariable{i1},'=value_iqc(',aloutvariable{i1},');'];
        exe_comment{sc}=str;
        sc=sc+1;
    end
end

sc=sc-1;


function [IQC_INFO]=add_blk(I,filename)

% --- check how many mudules to be added --- %
IQC_INFO=I;
name=[];


%--------------------------------------得到主層下的Inport/Outport數量
in_counter=size(find_system(filename,'SearchDepth',1,'BlockType','Inport'),1);
out_counter=size(find_system(filename,'SearchDepth',1,'BlockType','Outport'),1);


for nth=1:size(I,1),
    
    % --- extract iqc box's position --- %
    
    x1=I(nth,1);
    y1=I(nth,2);
    x2=I(nth,3);
    y2=I(nth,4);
    input_width=I(nth,6);
    output_width=I(nth,5);
    
    % --- 1. add inport block --- %
    if input_width~=0,
        IQC_INFO(nth,8)=in_counter+1;   % stor information into IQC_INFO
        
        c1=x2;
        c2=y2;
        c3=c1+20;
        c4=c2+20; % c1,c2,c3,c4 are block's coordinates %
        blk_pos=['[',num2str(c1),',',num2str(c2),',',num2str(c3),...
            ',',num2str(c4),']'];
        
        name=['MIn', num2str(in_counter+1)];
        add_block('built-in/Inport',[filename,'/',name],'Position',blk_pos,...
            'Port',num2str(in_counter+1),'PortWidth',num2str(input_width),...
            'SampleTime','-1','SignalType','auto','Interpolate','on');
        in_counter=in_counter+1;
    end
    
    % --- 2. add outport blocks --- %
    if output_width~=0,
        IQC_INFO(nth,7)=out_counter+1;
        
        c1=x1;
        c2=y2;
        c3=c1+20;
        c4=c2+20; % c1,c2,c3,c4 are block's coordinates %
        blk_pos=['[',num2str(c1),',',num2str(c2),',',num2str(c3),...
            ',',num2str(c4),']'];
        
        name=['MOut',num2str(out_counter+1)];
        add_block('built-in/Outport',[filename,'/',name],'Position',blk_pos,...
            'Port',num2str(out_counter+1),'PortWidth',num2str(output_width),...
            'SampleTime','-1','OutputWhenDisabled','held','InitialOutput','0');
        out_counter=out_counter+1;
    end
end

function linechang(LinearFileName,iqc_block,IQC,IQC_INFO)

line_item=find_system(LinearFileName,'FindAll','on','type','line');
line_src=get_param(line_item,'SrcBlock');
line_src_port=get_param(line_item,'SrcPort');
line_dst=get_param(line_item,'DstBlock');
line_dst_port=get_param(line_item,'DstPort');

delete_block(iqc_block)
% delete_block(L2block)


for n1=1:length(IQC.name)
    for n2=1:size(line_item,1)
        
        if strcmp(line_src{n2},IQC.name{n1}) && IQC_INFO(n1,8)~=0
            try
                delete_block(iqc_block{n1})
            end
            dstpoint=get_param(line_item(n2),'Points');
            in_point=get_param(find_system(LinearFileName,'Name',['MIn',num2str(IQC_INFO(n1,8))]),'Position');
            in_point_x=in_point{1}(1,3)+5;
            in_point_y=in_point{1}(1,4)-10;
            add_line(LinearFileName,[in_point_x,in_point_y;dstpoint(1,:)]);
            %             return;
        end
        %--------------------------------設定連接Outport
        if strcmp(line_dst{n2},IQC.name{n1}) && IQC_INFO(n1,7)~=0;
            out_point=get_param(find_system(LinearFileName,'Name',['MOut',num2str(IQC_INFO(n1,7))]),'Position');
            out_point_x=out_point{1}(1,1)-5;
            out_point_y=out_point{1}(1,2)+10;
            srcpoint=get_param(line_item(n2),'Points');
            try
                delete_block(iqc_block{n1})
            end
            add_line(LinearFileName,[srcpoint(size(srcpoint,1),:);out_point_x,out_point_y]);
            %             return;
        end
    end
end
% out=1;

function IQC_INFO=getportdim(I,filename)
IQC_INFO=I;
in_counter=1;
out_counter=1;

for i1=1:size(I,1)
    % outport signal pos
    if IQC_INFO(i1,5)~=0
        IQC_INFO(i1,9)=out_counter;
        IQC_INFO(i1,10)=out_counter+IQC_INFO(i1,5)-1;
        out_counter=out_counter+IQC_INFO(i1,5);
    end
    % inport signal pos
    if IQC_INFO(i1,6)~=0
        IQC_INFO(i1,11)=in_counter;
        IQC_INFO(i1,12)=in_counter+IQC_INFO(i1,6)-1;
        in_counter=in_counter+IQC_INFO(i1,6);
    end
end



function [A]=findpar1(line_buffer)

pos1=size(line_buffer,1);
pos2=size(line_buffer,2);
num=[];
A=[];
s=0;
for cur=pos1:pos2,
    if any('1234567890'==line_buffer(cur))
        num=[num,line_buffer(cur)];
        s=1;
    elseif s==1,
        s=0;
        A=[A,str2num(num)];
        num=[];
    end
end

function idx=findcell(m_cell,f_str)
nx=size(m_cell,2);
for i=1:nx
    if strcmp(m_cell{i},f_str)
        idx=i;
        return;
    end
end
idx=[];
return;