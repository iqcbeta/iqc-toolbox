function ext_sim_window(block)
setup(block);

function setup(block)

block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 11;

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions = -1;
block.InputPort(1).DirectFeedthrough = false;
block.InputPort(1).SamplingMode = 'sample';
block.InputPort(1).SampleTime = [-1 0];

block.OutputPort(1).Dimensions = -1;
block.OutputPort(1).SamplingMode = 'sample';
block.OutputPort(1).SampleTime = [0 0];

block.RegBlockMethod('SetInputPortDimensions',@SetInputPortDims);
block.RegBlockMethod('SetInputPortSampleTime',@SetInputSampleTime);
block.RegBlockMethod('SetOutputPortSampleTime',@SetOutputSampleTime);
block.RegBlockMethod('CheckParameters',@CheckParameters);
block.RegBlockMethod('Outputs', @Outputs);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Terminate', @Terminate);
block.RegBlockMethod('InitializeConditions',@InitializeConditions);

function InitializeConditions(block)
str=['text(5,85,''Sim IQC'');',...
    'plot(0,0,100,100,[84,94,4],[27,25,25],[94,84],[25,22],[15,15],',...
    '[5,97],[12,15],[86,97],[18,15],[86,97],[81,70,70,16,16,6],',...
    '[26,26,72,72,26,26],[71,71,15],[26,73,73])'];
block_path=gcb;
dot_pos=strfind(block_path,'/');
set_param(block_path(1:dot_pos(length(dot_pos))-1),'MaskDisplay',str);


function CheckParameters(block)
global ext_iqc_setvalues
block_path=gcb;
try
    if isempty(ext_iqc_setvalues.ext_sim_window.idx)
        i_pos=1;
    else
        i_pos=ext_iqc_setvalues.ext_sim_window.idx;
        findi=find_ipos(gcb);
        if isempty(findi)
            i_pos=i_pos+1;
        else
            i_pos=findi;
            return
        end
    end
catch err
    i_pos=1;
end
block_name=get_param(gcb,'name');
ext_iqc_setvalues.ext_sim_window.path{i_pos}=block_path;
ny=strfind(block_path,'/');
block_path=[block_path(1:ny-1),'_lin',block_path(ny:ny(end)-1)];
ext_iqc_setvalues.ext_sim_window.newpath{i_pos}=block_path;
ext_iqc_setvalues.ext_sim_window.name{i_pos}=block_name;
if ~isfield(ext_iqc_setvalues.ext_sim_window,'idx')
    ext_iqc_setvalues.ext_sim_window.idx=i_pos;
elseif i_pos>ext_iqc_setvalues.ext_sim_window.idx
    ext_iqc_setvalues.ext_sim_window.idx=i_pos;
end

% Simulink Parameters Check
max_del=block.DialogPrm(2).Data;
if ~strcmp(class(max_del), 'double')
    DAStudio.error('Maximum Delay must be numeric');
end
if max_del<0
    DAStudio.error('Maximum Delay must be positive number');
end
num_ran=block.DialogPrm(3).Data;
if ~strcmp(class(num_ran), 'double')
    DAStudio.error('Number of Random must be numeric');
end
if num_ran<0
    DAStudio.error('Number of Random must be positive number');
end
cus=block.DialogPrm(4).Data;
if ~strcmp(class(cus), 'double')
    DAStudio.error('Custom Delay must be numeric');
end
if any(cus<0)
    DAStudio.error('Custom Delay must be positive number');
end
if any(cus>max_del)
    DAStudio.error('Maximum Delay must be greater than Custom Delay');
end
var_nam=block.DialogPrm(5).Data;
if ~ischar(var_nam)
    DAStudio.error('Variable Name must be char');
end

% IQC Procedure Parameters Check
iqc1=block.DialogPrm(6).Data;
iqc2=block.DialogPrm(9).Data;
iqc3=block.DialogPrm(12).Data;
if ~iqc1 && ~iqc2 && ~iqc3
    DAStudio.error('The IQC procedure parameters set error');
    clear global ext_iqc_setvalues ABST
end
if iqc1
    iqc1_pole=block.DialogPrm(7).Data;
    if ~strcmp(class(iqc1_pole), 'double')
        DAStudio.error('iqc_delay Pole must be numeric');
    end
    iqc1_out=block.DialogPrm(8).Data;
    for i=1:length(iqc1_out)
        if ~isvarname(iqc1_out{i})
            DAStudio.error('iqc_delay Output Variable set error');
        end
    end
    setvalue(i_pos,'iqc_delay',iqc1_out,{'max_tim',max_del;...
        'pole',iqc1_pole});
end
if iqc2
    iqc2_pole=block.DialogPrm(10).Data;
    if ~strcmp(class(iqc2_pole), 'double')
        DAStudio.error('iqc_delay1 Pole must be numeric');
    end
    iqc2_out=block.DialogPrm(11).Data;
    for i=1:length(iqc2_out)
        if ~isvarname(iqc2_out{i})
            DAStudio.error('iqc_delay1 Output Variable set error');
        end
    end
    setvalue(i_pos,'iqc_delay1',iqc2_out,{'max_tim',max_del;...
        'pole',iqc2_pole});
end
if iqc3
    lis_rel=block.DialogPrm(13).Data;
    if ~ischar(lis_rel)
        DAStudio.error('custom IQC relationships must be char');
    end
    setvalue(i_pos,'iqc_custom',{'0'},{'lis_rel',lis_rel});
end

function Start(block)
global ext_iqc_setvalues
i_pos=find_ipos(gcb);
delay_time=block.DialogPrm(14).Data;
ext_iqc_setvalues.ext_sim_window.delay_time{i_pos}=delay_time;
inputdim=block.InputPort(1).Dimensions;
block_name=gcb;
dot_pos=strfind(gcb,'/');
n=length(dot_pos);
block_name=block_name(dot_pos(n-1)+1:dot_pos(n)-1);
fprintf('The "%s" Simulink Parameter\n',block_name)
if inputdim==size(delay_time,2)
    for i=1:size(delay_time,2)
        fprintf('The #%d Random Delay Time: %4.3f \n',...
            i,delay_time(1,i));
    end
elseif inputdim~=1 && size(delay_time,2)==1
    for i=1:inputdim
        fprintf('The #%d Random Delay Time: %4.3f \n',...
            i,delay_time(1,1));
    end
else
    error('The Number of Random ~= Input dimension')
end
var_nam=block.DialogPrm(5).Data;
if ~isempty(var_nam)
    eval([var_nam,'=delay_time;']);
    eval(['global ',var_nam,';']);
end

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
i_pos=find_ipos(gcb);
ext_iqc_setvalues.ext_sim_window.sampletime{i_pos}(1:2)=di(1);
if di(1)==0,
    block.InputPort(idx).Sampletime = [0 di(2)];
    block.OutputPort(idx).Sampletime = [0 0];
else
    block.InputPort(idx).Sampletime = di;
    block.OutputPort(idx).Sampletime = [di(1) 0];
end

function SetOutputSampleTime(block, idx, di)


function SetInputPortDims(block, idx, di)
global ext_iqc_setvalues
i_pos=find_ipos(gcb);
ext_iqc_setvalues.ext_sim_window.size{i_pos}(1:2)=di;
% Set compiled dimensions
block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = di;

% block_path=gcb;
% ny=strfind(block_path,'/');
% block_path=block_path(1:ny(end)-1);
% set_param(block_path,'num_ran',num2str(block.InputPort(1).Dimensions))


function Outputs(block)
i_pos=find_ipos(gcb);
global ext_iqc_setvalues
if ~isfield(ext_iqc_setvalues.ext_sim_window,'size')
    ext_iqc_setvalues.ext_sim_window.size{i_pos}(1:2)=...
        block.InputPort(1).Dimensions;
end
block.OutputPort(1).Data=zeros(block.OutputPort(1).Dimensions,1);

%endfunction

function Terminate(block)
iqctool('clear')
%end Terminate

function setvalue(i_pos,iqc_name,out_var,other)
global ext_iqc_setvalues
str=[];
for i=1:length(out_var)
    str=[str,out_var{i},','];
end
try
    ny1=size(ext_iqc_setvalues.ext_sim_window.iqctype{i_pos},2);
catch err
    ny1=0;
end
if ~strcmp(str(1:end-1),'0')
    ext_iqc_setvalues.ext_sim_window.iqctype{i_pos}{ny1+1}(1,:)={[iqc_name,...
        '.outvariable'],...
        str};
else
    ext_iqc_setvalues.ext_sim_window.iqctype{i_pos}{ny1+1}(1,:)={[iqc_name,...
        '.list_iqc_rel'],...
        other{1,2}};
    other=[];
end
for i=1:size(other,1)
    
    ny2=size(ext_iqc_setvalues.ext_sim_window.iqctype{i_pos}{ny1+1},1);
    ext_iqc_setvalues.ext_sim_window.iqctype{i_pos}{ny1+1}(ny2+1,:)=...
        {[iqc_name,'.',other{i,1}],other{i,2}};
end

function findi=find_ipos(gcb)
global ext_iqc_setvalues
i_pos=ext_iqc_setvalues.ext_sim_window.idx;
for i=1:i_pos,
    if strcmp(ext_iqc_setvalues.ext_sim_window.path(i),gcb)
        findi=i;
        return
    end
end
findi=[];