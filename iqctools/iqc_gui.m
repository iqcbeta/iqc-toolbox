function [gain,exe_comment]=iqc_gui(sys_name,datafile_name,type)
% function [gain,exe_comment]=iqc_gui(sys_name,datafile_name,type)
% This is the graohic interface function for the IQC ToolBox 
% 
% Created by C. Kao cykao@mit.edu  lasted modified: 10/27/1997
%                                  lasted modified: 06/06/1998 -- add some new functions
%                                  lasted modified: 07/09/1999 -- fix the SIMULINK3 compatability
%                                                                 problem


%
% corrected a problem with internal variable "comp", changed to "compstr"
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu


abst_init_iqc;   % initialize the abst environment 
global ABST
options=ABST.options;
fclose('all');

if nargin==0,
   disp('You have to build a SIMULINK block diagram')
   disp('using the blocks from "iqc_lib" before')
   disp('using this program. ')
   sys_name=input('Enter diagram "filename" now: ','s');
   disp(' ')
elseif nargin==1,
   datafile_name=[];
   isdatafile=0;
elseif nargin==2,
   isdatafile=1;
   type=1;
   dot_pos=(findstr(datafile_name,'.'));
   if ~isempty(dot_pos)
       datafile_name=datafile_name(1:dot_pos(1,1)-1);
   end
   datafile_name_full=[datafile_name,'.m'];
elseif nargin==3,
   isdatafile=1;
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

if isdatafile
   if ~exist(datafile_name_full)
      disp('The specified data file does not exist ...') 
      error('Or the file is not in the directory which is on the MATLAB search path')
   end
end

if options(5)==66 | isdatafile,
   keep_linfile=1;
else
   keep_linfile=0;
end

[exe_comment,ecctr]=gui_iqctbx(sys_name,keep_linfile);
exe_sys_name=[sys_name,'_iqcexe'];
exe_sys_name_full=[sys_name,'_iqcexe.m'];

fp=fopen(exe_sys_name_full,'wt');
fprintf(fp,['function gain=',exe_sys_name,'\n']);
fprintf(fp,['abst_init_iqc;\n']);

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
   eval(['abst_init_iqc;'])
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

% gain=eval(exe_sys_name);



function [exe_comment,ecctr,A,B,C,D,newfile,IQC]=gui_iqctbx(filename,keep_linfile)
% [exe_comment,ecctr,A,B,C,D,newfile,IQC]=gui_iqctbx(filename,keep_linfile)
%
% This function precess the SIMULink diagram, extract the
% information from iqc-blocks, replace them by Mux/Demux 
% blocks, and compute the A,B,C,D matrices of the resulting 
% linear part. Then it defines all signals and corresponding
% IQC relationships, and solves the L2 gain of the system by 
% using IQC Toolbox created by %
% Input parameters:
% -filename    file with graphic representation 
% of Simulink block diagram
%
% Output parameters:
% -gain        L2 gain of the system interested
% -A, B, C, D  system matrices of the linear part of the block diagram 
% -newfile     file with LTI blocks only
% 
% and also a structure IQC, which contains information of the iqc-relations
% of the signals in the block diagram
%
% structure of IQC:
% -IQC.name: cell contains names of all iqc-blocks
% -IQC.rela: nx2 structure contains information of iqc-blocks
%  first  column is the desctrption; i.e. relation 
%  of input and output of the block 
%  second column is the size of the block; 
%  ie, [n,m]: block has n outputs and m inputs

if nargin==1,
   keep_linfile=0;
end

global ABST
options=ABST.options;
global exe_comment
exe_comment={};           % exectuable comments cell array
ecctr=0;                  % exectuable comments cell array counter

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

% First round: Process the Simulink diagram
% 1. extract information contained in iqc-blocks
% 2. replace iqc-blocks by Mux and DeMux blocks
% 3. compute the A B C D matrices of the linear part
%    of the diagram by linmod.m function

if (options(5)==0)|(options(5)==777),
   disp(['extracting information from ', filename, ' ...'])
   disp(' ')
end

filenamemdl=[filename '.mdl'];
FID=fopen(filenamemdl,'rt');
if (FID==-1)
   fprintf('Can not open this file -- %s!\n',filename)
   return;
end

newfile=[filename '_lin'];
newfilemdl=[newfile '.mdl'];
if exist(newfilemdl),
   delete(newfilemdl);
end

Fid=fopen(newfilemdl,'wt');

% -- check if this is a simulink file -- %
line_temp=fgets(FID);
if ~strcmp(line_temp(1:7),'Model {'),
    disp('Error in using new_extr.m!');
    return;
else
    fprintf(Fid,line_temp); %write the first line into newfile%
end

% -- Main loop -- %

IQC_INFO=[];      % containing information of IQC block
have_done=0;      % index for adding Mux and DeMux mudules
iqc_counter=0;    % IQC block counter
block_counter=0;  % block counter

while ~feof(FID),
  line_temp=fgets(FID);

  % -- check file name and replace by new one -- %

  cor=findstr(line_temp,'%');
  if ~isempty(cor)
     count=0;
     strtmp=[line_temp(1:cor(1)-1),'%%'];
     for flcnt=2:length(cor)
         strtmp=[strtmp,line_temp(cor(flcnt-1)+1:cor(flcnt)-1),'%%'];
     end
     strtmp=[strtmp,line_temp(cor(length(cor))+1:length(line_temp))];
     line_temp=strtmp;
  end
    
  if ~isempty(findstr(line_temp,filename)),
      line_temp=strrep(line_temp,filename,newfile);
      fprintf(Fid,line_temp);
  

  % -- Find Block{ } and copy it to new file if it doesn't relate to IQC 
  %    If it is IQC block, extract the information contained in it  
  %    Also counter the number of Mux block and Demux block

  elseif ~isempty(findstr(line_temp,'Block {')),
     block_counter=block_counter+1;
     Buffer=[];
     Buffer=[Buffer,line_temp];
     iqc_index=0;
     while isempty(findstr(line_temp,'}')),
       line_temp=fgets(FID);
       Buffer=[Buffer,line_temp];
     % length=size(line_temp,2);

     % extracting the name of the block and check if the block is IQC block
   
       if size(line_temp,2)>11,
          if strcmp(line_temp(6:11),setstr([32 78 97 109 101 9])),
             cor=findstr(line_temp,'"');
             eval(['global used_name',num2str(block_counter)]);
             eval(['used_name',num2str(block_counter),...
                   '=line_temp(1,cor(1,1)+1:cor(1,2)-1);']); 
          end
          if ~isempty(findstr(line_temp,'S-Function')),
             iqc_index=1;
             iqc_counter=iqc_counter+1;
          end
       end
     end 
   
     if ~(iqc_index==1),       % Not a iqc block, simply copy 
         fprintf(Fid,Buffer);         
     else                      % Iqc block founded! extracting the information
         [I,name,fun]=iqc_info(Buffer);
         IQC_INFO=[IQC_INFO;I];
         IQC.name{iqc_counter}=name;
         IQC.rela{iqc_counter,1}=fun;
         IQC.rela{iqc_counter,2}=I(5:6);
     end

  % -- If Line{ } found, (it means all Block{ } have been exhausted), then
  %    add Inport and Outport into the file --  

   elseif ~isempty(findstr(line_temp,'Line {')), 
     if have_done==0, 
         IQC_INFO=add_blk(Fid,IQC_INFO,block_counter);
         have_done=1;
         fprintf(Fid,line_temp);
     else
        fprintf(Fid,line_temp);   
     end

  % -- Find Line{ } and copy it to new file if it doesn't relate to IQC
  %    and replace any iqc part by mux and demux -- 

  elseif ~isempty(findstr(line_temp,'SrcBlock')),
%     if ~isempty(findstr(line_temp,'iqc-block')),
        cor=findstr(line_temp,'"');
        name=line_temp(1,cor(1,1)+1:cor(1,2)-1);
        repi=0;
        for k=1:iqc_counter,
             if strcmp(IQC.name{k},name)
               rep=['In',num2str(IQC_INFO(k,7))];
               repi=1;
            end
        end
        if repi,
           line_temp=strrep(line_temp,name,rep);
        end
        fprintf(Fid,line_temp);

   elseif ~isempty(findstr(line_temp,'DstBlock')),
%      if ~isempty(findstr(line_temp,'iqc-block')),
        cor=findstr(line_temp,'"');
        name=line_temp(1,cor(1,1)+1:cor(1,2)-1);
        repi=0;
        for k=1:iqc_counter,
             if strcmp(IQC.name{k},name)
               rep=['Out',num2str(IQC_INFO(k,8))];
               repi=1;
            end
        end
        if repi,
           line_temp=strrep(line_temp,name,rep);
        end
        fprintf(Fid,line_temp);
%     else
%        fprintf(Fid,line_temp);
%     end
     
  % -- Nothing happened, simply copy it to new file -- %
  else 
     fprintf(Fid,line_temp);  
  end
end

fclose(FID);
fclose(Fid);


% Second Round: 
% Defind the signals and corresponding IQC relationships of
% the signals, and run the solver

if (options(5)==0)|(options(5)==777),
   disp('define the system and run solver ...')
end

exe_comment{ecctr+1}=['open_system(''',newfile,''');'];
exe_comment{ecctr+2}=['[A_iqcexe,B_iqcexe,C_iqcexe,D_iqcexe]=linmod2(''',newfile,''');'];
exe_comment{ecctr+3}=['close_system(''',newfile,''',1);'];
ecctr=ecctr+3;
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

exe_comment{ecctr+1}=['w_iqcexe=ss(A_iqcexe,B_iqcexe,C_iqcexe,D_iqcexe)*f_iqcexe;'];
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
ecctr=ecctr+1;


function [A,block_name,fun]=iqc_info(Buffer)
% function [A,block_name,fun]=iqc_info(Buffer)
% This function is used in graph_extr.m, for extracting the 
% information of iqc-blocks in the bolck diagram
% 
% -A : 1x6 matrix 
%      A(1:4) : position of the block
%      A(5:6) : size of the block; ie, A(5) is the number of the
%      output of the block, A(6) is the number of the input of the 
%      block  
% -block_name: the name of the block (a string)
% -fun: the description of the block (a string)

A=zeros(1,8);
Buffer_length=size(Buffer,2);
pos=1;
block_name=[];
fun=[];
blkname_complete=0;
fun_complete=0;
% blkname_gointo=0;
% fun_gointo=0;

while pos<=Buffer_length,
    line_temp=[];
    while ~(abs(Buffer(pos))==10),        % get a line from Buffer %
       line_temp=[line_temp,Buffer(pos)];
       pos=pos+1;
    end
    
    % --- parameter value takes more than one line ---
    if ~isempty(block_name) & blkname_complete==0
       if size(line_temp,2)>9 & all(abs(line_temp(1:9))==[9,9,9,32,32,32,32,32,32]),
          cor1=findstr(line_temp,'"')
          block_name=[block_name,line_temp(cor1(1,1)+1:cor1(1,2)-1)];       
       else 
          blkname_complete=1; 
       end
    end
    if ~isempty(fun) & fun_complete==0
       if size(line_temp,2)>9 & all(abs(line_temp(1:9))==[9,9,9,32,32,32,32,32,32]),
          cor1=findstr(line_temp,'"');
          fun=[fun,line_temp(cor1(1,1)+1:cor1(1,2)-1)];       
       else 
          fun_complete=1;
       end
    end

    % --- extracting block name --- %
    if size(line_temp,2)>11,
     if all(abs(line_temp(6:11))==[32 78 97 109 101 9]),
        cor=findstr(line_temp,'"');
        block_name=line_temp(1,cor(1,1)+1:cor(1,2)-1);
     end
    end

    % --- extracting block position --- %
    if ~isempty(findstr(line_temp,'Position'))
       cor1=findstr(line_temp,'[');
       cor2=findstr(line_temp,']');
       line_buffer=line_temp(1,cor1:cor2);
       A(1,1:4)=findpar1(line_buffer);
    end

    % --- extracting block order and related iqc function -- %
    if ~isempty(findstr(line_temp,'MaskValueString')),
       cor1=findstr(line_temp,'"');
       cor2=findstr(line_temp,'|');
       line_buffer=line_temp(1,cor1(1,1):cor2);
       order=findpar1(line_buffer);
       if size(order,2)==1,
          A(1,5)=order;
          A(1,6)=order;
       elseif size(order,2)==2,
          A(1,5:6)=order;
%       else
%          fprintf('Error occurs in extracing infromation from...
%                   %s',block_name);
%          break;
       end
       fun=line_temp(1,cor2(1,1)+1:cor1(1,2)-1);
    end
    
%    % --- parameter value takes more than one line ---
%    if ~isempty(block_name) & blkname_complete==0
%       if blkname_gointo==0
%          blkname_gointo=1;
%       else
%       if size(line_temp,2)>9 & all(abs(line_temp(1:9))==[9,9,9,32,32,32,32,32,32]),
%          cor1=findstr(line_temp,'"')
%          block_name=[block_name,line_temp(cor1(1,1)+1:cor1(1,2)-1)];       
%       else 
%          blkname_complete=1; 
%       end
%       end
%    end

%    if ~isempty(fun) & fun_complete==0
%       if fun_gointo==0
%          fun_gointo=1;
%       else
%       if size(line_temp,2)>9 & all(abs(line_temp(1:9))==[9,9,9,32,32,32,32,32,32]),
%          cor1=findstr(line_temp,'"');
%          fun=[fun,line_temp(cor1(1,1)+1:cor1(1,2)-1)];       
%       else 
%          fun_complete=1;
%       end
%       end
%    end
     pos=pos+1;
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


function [IQC_INFO]=add_blk(Fid,I,block_counter)

% --- check how many mudules to be added --- %
IQC_INFO=I;
n=size(I,1);;
in_counter=1;
out_counter=1;

for nth=1:n,

% --- extract iqc box's position --- %

x1=I(nth,1);
y1=I(nth,2);
x2=I(nth,3);
y2=I(nth,4);
input_width=I(nth,5);

% --- 1. add inport block --- %           
    name=['In',num2str(in_counter)];
    port=[num2str(in_counter)];
    for k=1:block_counter,                       
      eval(['global used_name',num2str(k)]);     % check if use the same name 
      eval(['compstr=used_name',num2str(k),';']);   % of existed blocks.
      if strcmp(compstr,name)
         in_counter=in_counter+1;
         name=['In',num2str(in_counter)];
         port=[num2str(in_counter)];
      end
    end
    IQC_INFO(nth,7)=in_counter;   % stor information into IQC_INFO    
    c1=x2;
    c2=y2;
    c3=c1+20;
    c4=c2+20; % c1,c2,c3,c4 are block's coordinates %
    blk_pos=['[',num2str(c1),',',num2str(c2),',',num2str(c3),...
             ',',num2str(c4),']']; 
    line_temp1=['      Name                   "',name,'"',setstr(10)];
    line_temp2=['      Position               ',blk_pos,setstr(10)];
    line_temp3=['      Port                   "',port,'"',setstr(10)];
    line_temp4=['      PortWidth              "',num2str(input_width),...
                '"',setstr(10)];
    fprintf(Fid,'    Block {\n');
    fprintf(Fid,'      BlockType              Inport\n');
    fprintf(Fid,line_temp1);
    fprintf(Fid,line_temp2);
    fprintf(Fid,line_temp3);
    fprintf(Fid,line_temp4);
    fprintf(Fid,'      SampleTime             "-1"\n');
    fprintf(Fid,'      DataType               "auto"\n');
    fprintf(Fid,'      SignalType             "auto"\n');
    fprintf(Fid,'      Interpolate            on\n');    
    fprintf(Fid,'    }\n');


 % --- 2. add outport blocks --- %
    name=['Out',num2str(out_counter)];
    port=[num2str(out_counter)];
    for k=1:block_counter,                      
      eval(['global used_name',num2str(k)]);     % check if use the same name 
      eval(['compstr=used_name',num2str(k),';']);   % of existed blocks.
      if strcmp(compstr,name)
         out_counter=out_counter+1;
         name=['Out',num2str(out_counter)];
         port=[num2str(out_counter)];
      end
    end
    IQC_INFO(nth,8)=out_counter;
    c1=x1;
    c2=y2;
    c3=c1+20;
    c4=c2+20; % c1,c2,c3,c4 are block's coordinates %
    blk_pos=['[',num2str(c1),',',num2str(c2),',',num2str(c3),...
             ',',num2str(c4),']']; 
    line_temp1=['      Name                   "',name,'"',setstr(10)];
    line_temp2=['      Position               ',blk_pos,setstr(10)];
    line_temp3=['      Port                   "',port,'"',setstr(10)];   
    fprintf(Fid,'    Block {\n');
    fprintf(Fid,'      BlockType              Outport\n');
    fprintf(Fid,line_temp1);
    fprintf(Fid,line_temp2);
    fprintf(Fid,line_temp3);
    fprintf(Fid,'      OutputWhenDisabled     held\n');
    fprintf(Fid,'      InitialOutput          "0"\n');
    fprintf(Fid,'    }\n');

  in_counter=in_counter+1;
  out_counter=out_counter+1;
end










