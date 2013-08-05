function [gain,exe_comment]=iqc_gui(sys_name,datafile_name,type,sol)
%function [gain,exe_comment]=iqc_gui(sys_name,datafile_name,type,sol)
% if:       sys_name = 'newtest1'
%           sol = 'sdpt3' use solver: sdpt3 and call YALMIP
%
% command >> iqc_gui('newtest1',[],[],'sdpt3')
%
%



% options=[0 200 0 0 777];
options=[0 200 0 0 66];
fclose('all');

%%
% 判斷輸入程序
if nargin==0,
    disp('You have to build a SIMULINK block diagram')
    disp('using the blocks from "iqc_lib" before')
    disp('using this program. ')
    sys_name=input('Enter diagram "filename" now: ','s');
    disp(' ')
elseif nargin==1,
    datafile_name=[];
    isdatafile=0;
    sol=[];
elseif nargin==2
    isdatafile=1;
    type=1;
    sol=[];
    datafile=datafile_name;
    dot_pos=(findstr(datafile_name,'.'));
    if ~isempty(dot_pos)
        datafile_name=datafile_name(1:dot_pos(1,1)-1);
    end
    datafile_name_full=[datafile_name,'.m'];
elseif nargin==3
    isdatafile=1;
    datafile=datafile_name;
    sol=[];
    dot_pos=(findstr(datafile_name,'.'));
    if ~isempty(dot_pos)
        datafile_name=datafile_name(1:dot_pos(1,1)-1);
    end
    if type==1,
        datafile_name_full=[datafile_name,'.m'];
    elseif type==2,
        datafile_name_full=[datafile_name,'.mat'];
    else
        error('Data file can be only MATLAB M-file or Mat-datafile')
    end
elseif nargin==4 && isempty(datafile_name) && isempty(type)
    datafile_name=[];
    isdatafile=0;
elseif nargin==4
    isdatafile=1;
    datafile=datafile_name;
    dot_pos=(findstr(datafile_name,'.'));
    if ~isempty(dot_pos)
        datafile_name=datafile_name(1:dot_pos(1,1)-1);
    end
    if type==1,
        datafile_name_full=[datafile_name,'.m'];
    elseif type==2,
        datafile_name_full=[datafile_name,'.mat'];
    else
        error('Data file can be only MATLAB M-file or Mat-datafile')
    end
else
    error('Too many inputs')
end


%%
if isdatafile
    if ~exist(datafile_name_full)
        disp('The specified data file does not exist ...')
        error('Or the file is not in the directory which is on the MATLAB search path')
    end
end

if options(5)==777 || isdatafile,
    keep_linfile=1;
else
    keep_linfile=0;
end

[exe_comment,ecctr]=gui_iqctbx(sys_name,keep_linfile,sol);

exe_sys_name=[sys_name,'_iqcexe'];
exe_sys_name_full=[sys_name,'_iqcexe.m'];

fp=fopen(exe_sys_name_full,'wt');
fprintf(fp,['function gain=',exe_sys_name,'\n']);
% fprintf(fp,['abst_init_iqc;\n']);


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
    %     eval(['abst_init_iqc;'])
    for forlp_counter=1:ecctr
        eval(exe_comment{forlp_counter})
    end
    if options(5)<=100,
        delete(exe_sys_name_full);
    end
else
    disp(' ')
    strtmp1='iqc_gui.m has processed the graphic file ';
    strtmp2='... and produced the corresponding executable M-file ';
    disp([strtmp1,sys_name,'.mdl'])
    disp([strtmp2,exe_sys_name_full])
end

%---------------------2013/04/15 delete old file ex: *.mdl.r*
old_sim_name={'r14','r14sp1','r14sp2','r14sp3','r2006a','r2006b',...
    'r2007a','r2007b','r2008a','r2008b','r2009a','r2009b',...
    'r2010a','r2010b','r2011a','r2011b','r2012a','r2012b','unknown_version'};
for i1=1:length(old_sim_name)
    if exist([sys_name,'_lin.mdl.',old_sim_name{i1}])
        delete([sys_name,'_lin.mdl.',old_sim_name{i1}]);
    end
end

%%
function [exe_comment,ecctr,A,B,C,D,newfile,IQC]=gui_iqctbx(filename,keep_linfile,sol)
%%
if nargin==1,
    keep_linfile=0;
end
global ABST
clear ABST
options=[0 200 0 0 0];

global exe_comment
exe_comment={};           % executable comments cell array
ecctr=0;                  % executable comments cell array counter

if (options(5)==0)|(options(5)==777),
    disp(' ')
    disp(' ')
    disp('*************************************************')
    disp(['iqc_gui: processing diagram "' filename '.mdl" ...'])
    disp('*************************************************')
    disp(' ')
end

global IQC
IQC.name={};
IQC.rela={};

if (options(5)==0)|(options(5)==777),
    disp(['extracting information from ', filename, ' ...'])
    disp(' ')
end
oldfilemdl=[filename '.mdl'];
newfile=[filename '_lin'];
newfilemdl=[newfile '.mdl'];

file_path=which(oldfilemdl);
dot_pos=(findstr(file_path,oldfilemdl));
cd(file_path(1:dot_pos-1));

if exist(newfilemdl),
    delete(newfilemdl);
end


copyfile(oldfilemdl,newfilemdl)

load_system(newfile);

IQC_INFO=[];      % containing information of IQC block
have_done=0;      % index for adding Mux and DeMux mudules
iqc_counter=0;    % IQC block counter
block_counter=0;  % block counter

%-------------------------------找IQC block 資訊
A=zeros(1,8);
fun=[];
%------------------------------block_counter需扣除subsystem個數及原系統
All_block=find_system(newfile);
Sub_block=find_system(newfile,'BlockType','SubSystem');
block_counter=size(All_block,1)-size(Sub_block,1)-1;

% -------------------------------------判斷為連續系統or離散系統
dis_block=find_system(newfile,'BlockType','DiscreteTransferFcn');
dis_block_sampletime=get_param(dis_block,'SampleTime');
dis_block=find_system(newfile,'BlockType','DiscreteZeroPole');
dis_block_sampletime=[dis_block_sampletime;get_param(dis_block,'SampleTime')];
dis_block=[dis_block;find_system(newfile,'BlockType','DiscreteStateSpace')];
dis_block_sampletime=[dis_block_sampletime;get_param(dis_block,'SampleTime')];
for i=1:(size(dis_block_sampletime,1)-1)
    if ~strcmp(dis_block_sampletime{i,1},dis_block_sampletime{i+1,1})
        error('sampletime error')
    end
end
if isempty(dis_block_sampletime)
    signal_type='continuous';
else
    signal_type='discrete';
    signal_sampletime=dis_block_sampletime{1,1};
end

% ------------------------------------找IQC方塊
iqc_block=find_system(newfile,'FunctionName','IQC');
iqc_block_name=get_param(iqc_block,'Name');
iqc_block_position=get_param(iqc_block,'Position');
iqc_block_MaskValueString=get_param(iqc_block,'MaskValueString');
iqc_counter=size(iqc_block,1);

for n=1:iqc_counter
    %--------------------------------記錄參數
    line_temp=iqc_block_MaskValueString{n};
    cor1=findstr(line_temp,'|');
    line_buffer=line_temp(1,1:cor1);
    order=findpar1(line_buffer);
    if size(order,2)==1,
        A(1,5)=order;
        A(1,6)=order;
    elseif size(order,2)==2,
        A(1,5:6)=order;
    end
    
    %--------------------------------記錄位置
    A(1,1:4)=iqc_block_position{n};
    IQC_INFO=[IQC_INFO;A];
    fun=line_temp(1,cor1(1,1)+1:size(line_temp,2));
    
    IQC.name{n}=iqc_block_name{n};
    IQC.rela{n,1}=fun;
    IQC.rela{n,2}=A(5:6);
    
end

%% 11/2
% delete main floor inport and outport
main_in=find_system(newfile,'SearchDepth',1,'BlockType','Inport');
delete_block(main_in);
main_out=find_system(newfile,'SearchDepth',1,'BlockType','Outport');
delete_block(main_out);
save_system(newfile,[],'OverwriteIfChangedOnDisk',true,'BreakAllLinks',true);

% additive Inport and Outport
IQC_INFO=add_blk(IQC_INFO,newfile);
save_system(newfile,[],'OverwriteIfChangedOnDisk',true);


% change the connect line
linechang(newfile,iqc_block,IQC,IQC_INFO);
save_system(newfile,[],'OverwriteIfChangedOnDisk',true);

close_system(newfile);

% IQC_INFO=getportdim(IQC_INFO,newfile);

% if length(unique(sampletimeset))~=1
%     error('the block sample time does not match')
% end

%%
% %--------------------------------增加輸入輸出
% IQC_INFO=add_blk(IQC_INFO,block_counter,newfile,All_block,Sub_block,iqc_block_name);

% %--------------------------------更改連接線
% line_item=find_system(newfile,'FindAll','on','type','line');
% line_src=get_param(line_item,'SrcBlock');
% line_dst=get_param(line_item,'DstBlock');
% inn=0;
% connet_line=[];
% for n1=1:size(iqc_block_name,1)
%     for n2=1:size(line_item,1)
%         %--------------------------------設定連接Inport
%         if strcmp(line_src{n2},iqc_block_name{n1});
%             Dstname=get_param(line_item(n2),'DstBlock');
%             dstpoint=get_param(line_item(n2),'Points');
%             in_point=get_param(find_system(newfile,'Name',['MIn',num2str(IQC_INFO(n1,8))]),'Position');
%             in_point_x=in_point{1}(1,3)+5;
%             in_point_y=in_point{1}(1,4)-10;
%             %----------------------------判斷分歧
%             if ~isempty(Dstname)
%                 delete_line(line_item(n2));
%                 add_line(newfile,[in_point_x,in_point_y;dstpoint(2:size(dstpoint,1),:)]);
%                 inn=inn+1;
%             else
%                 dstpoint=get_param(line_item(n2),'Points');
%                 %------------------------不能直接刪，若直接刪則分歧線也會刪除
%                 %------------------------先判斷分歧線連到哪
%                 for n3=1:size(line_item,1)
%                     l=get_param(line_item(n3),'Points');
%                     if dstpoint(size(dstpoint,1),:)==l(1,:)
%                         inn=inn+1;
%                         connet_line{inn}=get_param(line_item(n3),'Points');
%                     end
%                 end
%                 delete_line(line_item(n2));
%                 add_line(newfile,[in_point_x,in_point_y;dstpoint(2:size(dstpoint,1),:)]);
%                 for n4=1:inn
%                     add_line(newfile,[connet_line{n4}]);
%                 end
%                 inn=1;
%             end
%         end
%         %--------------------------------設定連接Outport
%         if strcmp(line_dst{n2},iqc_block_name{n1});
%             srcpoint=get_param(line_item(n2),'Points');
%             out_point=get_param(find_system(newfile,'Name',['MOut',num2str(IQC_INFO(n1,8))]),'Position');
%             out_point_x=out_point{1}(1,1)-5;
%             out_point_y=out_point{1}(1,2)+10;
%             delete_line(line_item(n2));
%             add_line(newfile,[srcpoint(1:size(srcpoint,1)-1,:);out_point_x,out_point_y]);
%         end
%     end
% end

% delete_block(iqc_block);
% save_system(newfile);
% close_system(newfile);

% Second Round:
% Defind the signals and corresponding IQC relationships of
% the signals, and run the solver
if (options(5)==0)|(options(5)==777),
    disp('define the system and run solver ...')
end
ecctr=1;

if strcmp(signal_type,'continuous')
    exe_comment{ecctr}=['abst_init_iqc;'];
else
    exe_comment{ecctr}=['abst_init_iqc(',signal_sampletime,');'];
end
ecctr=ecctr+1;

if ~isempty(sol)
    exe_comment{ecctr}=['setlmioptions(''',sol,''');'];
    ecctr=ecctr+1;
end


exe_comment{ecctr}=['open_system(''',newfile,''');'];
if strcmp(signal_type,'continuous');
    exe_comment{ecctr+1}=['[A_iqcexe,B_iqcexe,C_iqcexe,D_iqcexe]=linmod2(''',newfile,''');'];
else
    exe_comment{ecctr+1}=['[A_iqcexe,B_iqcexe,C_iqcexe,D_iqcexe]=dlinmod(''',newfile,''',',signal_sampletime,');'];
end
exe_comment{ecctr+2}=['close_system(''',newfile,''',1);'];
ecctr=ecctr+2;

if ~keep_linfile & options(5)<=100
    
    exe_comment{ecctr+1}=['delete(''',newfilemdl,''');'];
    ecctr=ecctr+1;
    
end

% --- count how many input does the problem have ---
fsize=0;
for forlp_counter=1:size(IQC.rela,1)
    fsize=fsize+IQC.rela{forlp_counter,2}(1,1);
end

% --- Find out the inputs which we are about to analysis and define
%     properties of inputs

f_counter=1;
w_counter=1;
special_input=0;

for forlp_counter=1:size(IQC.rela,1)
    COR1=findstr(IQC.rela{forlp_counter,1},'IN=');
    COR2=findstr(IQC.rela{forlp_counter,1},'OUT=');
    if ~isempty(COR1) & ~isempty(COR2)
        % get the index of components of the signal f which relate
        % to this iqc-block
        performance_block=forlp_counter; % mark the performance block
        fc1=num2str(f_counter);
        fc2=num2str(f_counter+IQC.rela{forlp_counter,2}(1)-1);
        OUT=['f_iqcexe(',fc1,':',fc2,')'];
        COR3=findstr(IQC.rela{forlp_counter,1},'|');
        if ~isempty(COR3)
            tempstr1=IQC.rela{forlp_counter,1}(COR3+1:size(IQC.rela{forlp_counter,1},2));
            if ~isempty(tempstr1)
                special_input=1;              % non-standard performance block founded!
                COR4=findstr(tempstr1,'=');
                if isempty(COR4)
                    error('Syntax error inside the performance box')
                elseif size(COR4,2)~=1
                    error('More than one properties specified in the performance box')
                else
                    tempstr2=tempstr1(1:(COR4-1));
                    COR5=findstr(tempstr2,',');
                    len=size(tempstr2,2);
                    if ~isempty(COR5) & strcmp(tempstr2(1,1),'[') & strcmp(tempstr2(1,len),']')
                        % extra outputs required
                        for forlp_counter=1:size(COR5,2)
                            if forlp_counter==size(COR5,2)
                                tempstr4=tempstr2(COR5(forlp_counter)+1:size(tempstr2,2)-1);
                                if isempty(tempstr4) | ~isempty(findstr(tempstr4,' '))
                                    error('Syntax error in performance block')
                                else
                                    exe_comment{ecctr+1}=['global ',tempstr4];
                                    ecctr=ecctr+1;
                                end
                            else
                                tempstr4=tempstr2(COR5(forlp_counter)+1:COR5(forlp_counter+1)-1);
                                if isempty(tempstr4) | ~isempty(findstr(tempstr4,' '))
                                    error('Syntax error in performance block')
                                else
                                    exe_comment{ecctr+1}=['global ',tempstr4];
                                    ecctr=ecctr+1;
                                end
                            end
                        end
                        tempstr3=['[f_iqcexe0',tempstr2(COR5(1,1):size(tempstr2,2))];
                    else
                        tempstr3='f_iqcexe0';
                    end
                    exe_comment{ecctr+1}=[tempstr3,tempstr1(COR4:size(tempstr1,2))];
                    ecctr=ecctr+1;
                end
            end
        end
        
        % get the index of components of the signal w which relates
        % to this iqc-block
        wc1=num2str(w_counter);
        wc2=num2str(w_counter+IQC.rela{forlp_counter,2}(2)-1);
        IN=['w_iqcexe(',wc1,':',wc2,')'];
    end
    f_counter=f_counter+IQC.rela{forlp_counter,2}(1);
    w_counter=w_counter+IQC.rela{forlp_counter,2}(2);
end

% define the input and output signals

exe_comment{ecctr+1}=['fsize=size(B_iqcexe,2);'];
ecctr=ecctr+1;

if special_input==1  % some inputs have specicified properties (have been defined)
    exe_comment{ecctr+1}=['fc1=',fc1,';','fc2=',fc2,';'];
    ecctr=ecctr+1;
    fc1=str2num(fc1);
    fc2=str2num(fc2);
    if fc1==1 & fc2==fsize, % the defined inputs are total inputs of this system
        exe_comment{ecctr+1}=['f_iqcexe=f_iqcexe0;'];
        ecctr=ecctr+1;
    elseif fc1==1 & fc2~=fsize,
        exe_comment{ecctr+1}=['f_iqcexe1=signal(fsize-fc2);'];
        exe_comment{ecctr+2}=['f_iqcexe=[f_iqcexe0;f_iqcexe1];'];
        ecctr=ecctr+2;
    elseif fc1~=1 & fc2==fsize,
        exe_comment{ecctr+1}=['f_iqcexe1=signal(fc1-1);'];
        exe_comment{ecctr+2}=['f_iqcexe=[f_iqcexe1;f_iqcexe0];'];
        ecctr=ecctr+2;
    else
        exe_comment{ecctr+1}=['f_iqcexe1=signal(fc1-1);'];
        exe_comment{ecctr+2}=['f_iqcexe2=signal(fsize-fc2);'];
        exe_comment{ecctr+3}=['f_iqcexe=[f_iqcexe1;f_iqcexe0;f_iqcexe2];'];
        ecctr=ecctr+3;
    end
else
    exe_comment{ecctr+1}=['f_iqcexe=signal(fsize);'];  % no inputs was defined before
    ecctr=ecctr+1;
end

% LTI system outputs, whose components are inputs of iqc-blocks
if strcmp(signal_type,'continuous');
    exe_comment{ecctr+1}=['w_iqcexe=ss(A_iqcexe,B_iqcexe,C_iqcexe,D_iqcexe)*f_iqcexe;'];
else
    exe_comment{ecctr+1}=['w_iqcexe=ss(A_iqcexe,B_iqcexe,C_iqcexe,D_iqcexe,',signal_sampletime,')*f_iqcexe;'];
end
ecctr=ecctr+1;

f_counter=1;         % reset counter
w_counter=1;         % reset counter

% --- Define IQC relationship ---

for forlp_counter=1:size(IQC.rela,1)
    if forlp_counter~=performance_block
        COR1=findstr(IQC.rela{forlp_counter,1},'=');
        COR1_n=size(COR1,2);
        tempstr1=IQC.rela{forlp_counter,1}(1:(COR1(1,1)-1));
        COR4=findstr(tempstr1,',');
        len=size(tempstr1,2);
        if ~isempty(COR4) & strcmp(tempstr1(1,1),'[') & strcmp(tempstr1(1,len),']')
            % extra outputs required
            tempstr3=tempstr1(2:COR4(1,1)-1);
            if isempty(tempstr3) | ~isempty(findstr(tempstr3,' '))
                error(['Syntax error in ',IQC.name{forlp_counter},' block'])
            end
            for forlp_counter1=1:size(COR4,2)
                if forlp_counter1==size(COR4,2)
                    tempstr2=tempstr1(COR4(forlp_counter1)+1:size(tempstr1,2)-1);
                    if isempty(tempstr2) | ~isempty(findstr(tempstr2,' '))
                        error(['Syntax error in ',IQC.name{forlp_counter},' block'])
                    else
                        exe_comment{ecctr+1}=['global ',tempstr2];
                        ecctr=ecctr+1;
                    end
                else
                    tempstr2=tempstr1(COR4(forlp_counter1)+1:COR4(forlp_counter1+1)-1);
                    if isempty(tempstr2) | ~isempty(findstr(tempstr2,' '))
                        error(['Syntax error in ',IQC.name{forlp_counter},' block'])
                    else
                        exe_comment{ecctr+1}=['global ',tempstr2];
                        ecctr=ecctr+1;
                    end
                end
            end
            extra_output=1;
        else
            extra_output=0;
        end
        
        len=size(IQC.rela{forlp_counter,1},2);
        tempstr4=IQC.rela{forlp_counter,1}(COR1(1,COR1_n)+1:len);
        len=size(tempstr4,2);
        COR2=findstr(tempstr4,'(');
        COR3=findstr(tempstr4,',');
        if isempty(COR3) % there is only one argument in iqc relation
            % description; it MUST be input of the block.
            % REPLACE it by the correct siginal defined here
            
            % get the index of components of the signal f which relate
            % to this iqc-block
            fc1=num2str(f_counter);
            fc2=num2str(f_counter+IQC.rela{forlp_counter,2}(1)-1);
            
            % get the index of components of the signal w which relate
            % to this iqc-block
            wc1=num2str(w_counter);
            wc2=num2str(w_counter+IQC.rela{forlp_counter,2}(2)-1);
            
            % construct the correct statement and evalute it!
            tp1=['f_iqcexe(',fc1,':',fc2,')'];
            tp2=['w_iqcexe(',wc1,':',wc2,')'];
            if extra_output==1
                exe_comment{ecctr+1}=['ub_',num2str(forlp_counter),'=',tp1,';'];
                exe_comment{ecctr+2}=[tempstr1,'=',tempstr4(1:COR2),tp2,');'];
                exe_comment{ecctr+3}=['ub_',num2str(forlp_counter),'==',tempstr3,';'];
                ecctr=ecctr+3;
            else
                exe_comment{ecctr+1}=[tp1,'==',tempstr4(1:COR2),tp2,');'];
                ecctr=ecctr+1;
            end
            
        else             % there is more than one arguments in iqc
            % relation; replace the first one by the
            % correct signal defined here.
            
            % get the index of components of the signal f which relate
            % to this iqc-block
            fc1=num2str(f_counter);
            fc2=num2str(f_counter+IQC.rela{forlp_counter,2}(1)-1);
            
            % get the index of components of the signal w which relate
            % to this iqc-block
            wc1=num2str(w_counter);
            wc2=num2str(w_counter+IQC.rela{forlp_counter,2}(2)-1);
            
            % construct the correct statement and evalute it!
            
            tp1=['f_iqcexe(',fc1,':',fc2,')'];
            tp2=['w_iqcexe(',wc1,':',wc2,')'];
            tp3=tempstr4(COR3:len);
            if extra_output==1
                exe_comment{ecctr+1}=['ub_',num2str(forlp_counter),'=',tp1,';'];
                exe_comment{ecctr+2}=[tempstr1,'=',tempstr4(1:COR2),tp2,tp3];
                exe_comment{ecctr+3}=['ub_',num2str(forlp_counter),'==',tempstr3,';'];
                ecctr=ecctr+3;
            else
                exe_comment{ecctr+1}=[tp1,'==',tempstr4(1:COR2),tp2,tp3];
                ecctr=ecctr+1;
            end
        end
    end
    f_counter=f_counter+IQC.rela{forlp_counter,2}(1);
    w_counter=w_counter+IQC.rela{forlp_counter,2}(2);
end

% Define the final 'run solver' comment
exe_comment{ecctr+1}=['gain=iqc_gain_tbx(',OUT,',',IN,');'];
ecctr=ecctr+2;

ecctr=ecctr-1;





% function [IQC_INFO]=add_blk(I,block_counter,newfile,All_block,Sub_block,iqc_block_name)
%
% % --- check how many mudules to be added --- %
% IQC_INFO=I;
% name=[];
% Sub_block_name=get_param(Sub_block,'Name');
% % -------------------------------------判斷IQC方塊所在位置是否在子方塊內
% % -------------------------------------如果有則題示使用者需將IQC拉出
% for i1=1:size(All_block,1)
%     for i2=1:size(Sub_block_name,1)
%         for i3=1:size(iqc_block_name,1)
%             if findstr(All_block{i1},[Sub_block_name{i2},'/',iqc_block_name{i3}])
%                 error('The IQC block must be in the main floor...')
%                 % error('請將IQC方塊拉出到主層')
%             end
%         end
%     end
% end
%
%
% %--------------------------------------得到主層下的Inport/Outport數量
% in_counter=size(find_system(newfile,'BlockType','Inport'),1);
% out_counter=size(find_system(newfile,'BlockType','Outport'),1);
%
% for i1=1:size(Sub_block_name,1)
%     sub_port=get_param(Sub_block,'Ports');
%
%     in_counter=in_counter-sub_port{1}(1);
%     out_counter=out_counter-sub_port{1}(2);
% end
%
% for nth=1:size(I,1),
%
%     % --- extract iqc box's position --- %
%
%     x1=I(nth,1);
%     y1=I(nth,2);
%     x2=I(nth,3);
%     y2=I(nth,4);
%     input_width=I(nth,5);
%     output_width=I(nth,6);
%
%     % --- 1. add inport block --- %
%     IQC_INFO(nth,7)=in_counter+1;   % stor information into IQC_INFO
%     name=['MIn', num2str(in_counter+1)];
%     c1=x2;
%     c2=y2;
%     c3=c1+20;
%     c4=c2+20; % c1,c2,c3,c4 are block's coordinates %
%     blk_pos=['[',num2str(c1),',',num2str(c2),',',num2str(c3),...
%         ',',num2str(c4),']'];
%     add_block('built-in/Inport',[newfile,'/',name],'Position',blk_pos,...
%         'Port',num2str(in_counter+1),'PortWidth',num2str(input_width),...
%         'SampleTime','-1','SignalType','auto','Interpolate','on');
%
%     % --- 2. add outport blocks --- %
%     IQC_INFO(nth,8)=out_counter+1;
%     name=['MOut',num2str(out_counter+1)];
%     c1=x1;
%     c2=y2;
%     c3=c1+20;
%     c4=c2+20; % c1,c2,c3,c4 are block's coordinates %
%     blk_pos=['[',num2str(c1),',',num2str(c2),',',num2str(c3),...
%         ',',num2str(c4),']'];
%     add_block('built-in/Outport',[newfile,'/',name],'Position',blk_pos,...
%         'Port',num2str(out_counter+1),'PortWidth',num2str(output_width),...
%         'SampleTime','-1','OutputWhenDisabled','held','InitialOutput','0');
%
%
%     in_counter=in_counter+1;
%     out_counter=out_counter+1;
% end


%% 11/2
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
    input_width=I(nth,5);
    output_width=I(nth,6);
    
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

%%
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

