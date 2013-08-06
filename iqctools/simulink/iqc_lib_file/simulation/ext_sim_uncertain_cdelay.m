function ext_sim_uncertain_cdelay(block)
setup(block);

function setup(block)

block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 10;

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions = -1;
block.InputPort(1).DirectFeedthrough = false;
block.InputPort(1).SamplingMode = 'sample';
block.InputPort(1).SampleTime = [-1 0];

block.OutputPort(1).Dimensions = -1;
block.OutputPort(1).SamplingMode = 'sample';
block.OutputPort(1).SampleTime = [-1 0];

block.RegBlockMethod('SetInputPortDimensions',@SetInputPortDims);
block.RegBlockMethod('SetInputPortSampleTime',@SetInputSampleTime);
block.RegBlockMethod('SetOutputPortSampleTime',@SetOutputSampleTime);
block.RegBlockMethod('CheckParameters',@CheckParameters);
block.RegBlockMethod('Outputs', @Outputs);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Terminate', @Terminate);
block.RegBlockMethod('InitializeConditions',@InitializeConditions);

function InitializeConditions(block)
str=['disp(''\fontsize{14}exp(-T*s)-1'', ''texmode'',''on'');',...
    'port_label(''input'',1,''in'');',...
    'port_label(''output'',1,''out'');',...
    'port_label(''output'',2,''T'');'];

block_path=gcb;
dot_pos=strfind(block_path,'/');
set_param(block_path(1:dot_pos(length(dot_pos))-1),'MaskDisplay',str);

function CheckParameters(block)
global ext_iqc_setvalues
% block_path=gcb;
% try
%     if isempty(ext_iqc_setvalues.ext_sim_uncertain_cdelay.idx)
%         pos=1;
%     else
%         pos=ext_iqc_setvalues.ext_sim_uncertain_cdelay.idx;
%         findi=find_ipos(gcb);
%         if isempty(findi)
%             pos=pos+1;
%         else
%             pos=findi;
%         end
%     end
% catch err
%     pos=1;
% end
% block_name=get_param(gcb,'name');
% ext_iqc_setvalues.ext_sim_uncertain_cdelay.path{pos}=block_path;
% ny=strfind(block_path,'/');
% block_path=[block_path(1:ny-1),'_lin',block_path(ny:ny(end)-1)];
% ext_iqc_setvalues.ext_sim_uncertain_cdelay.newpath{pos}=block_path;
% ext_iqc_setvalues.ext_sim_uncertain_cdelay.name{pos}=block_name;
% if ~isfield(ext_iqc_setvalues.ext_sim_uncertain_cdelay,'idx')
%     ext_iqc_setvalues.ext_sim_uncertain_cdelay.idx=pos;
% elseif pos>ext_iqc_setvalues.ext_sim_uncertain_cdelay.idx
%     ext_iqc_setvalues.ext_sim_uncertain_cdelay.idx=pos;
% end

% Simulink Parameters Check
max_del=block.DialogPrm(2).Data;
if ~isa(max_del,'double')
    disp_str(9,'maximum delay','numeric')
end
if max_del<0
    disp_str(9,'maximum delay','positive number')
end
num_ran=block.DialogPrm(3).Data;
if ~isa(num_ran,'double')
    disp_str(9,'number of random','numeric')
end
if num_ran<0
    disp_str(9,'number of random','positive number')
end
cus=block.DialogPrm(4).Data;
if ~isa(cus,'double')
    disp_str(9,'custom delay','numeric')
end
if any(cus<0)
    disp_str(9,'custom delay','positive number')
end
if any(cus>max_del)
    disp_str(39,'"maximum delay"','<','"custom delay"')
end

% IQC Procedure Parameters Check
block_path=gcb;
ny=strfind(block_path,'/');
block_newpath=[block_path(1:ny-1),'_lin',block_path(ny:ny(end)-1)];
block_name=get_param(gcb,'name');

glparam_name='ext_sim_uncertain_cdelay';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','path','newpath'};
alparam={block_name,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);

ext_iqc_setvalues.ext_sim_uncertain_cdelay.iqctype{pos}='';

iqc1=block.DialogPrm(5).Data;
iqc2=block.DialogPrm(8).Data;
if ~iqc1 && ~iqc2
    DAStudio.error('Simulink:block:invalidParameter');
    iqctool('simclear')
end
if iqc1
    iqc1_pole=block.DialogPrm(6).Data;
    if ~isa(iqc1_pole,'double')
        disp_str(9,'"iqc_cdelay pole"','numeric')
    end
    iqc1_out=block.DialogPrm(7).Data;
    for i=1:length(iqc1_out)
        if ~isvarname(iqc1_out{i})
            DAStudio.error('Simulink:block:invalidParameter');
        end
    end
    setvalue(pos,'iqc_cdelay',iqc1_out,{'max_tim',max_del;...
        'pole',iqc1_pole});
end
if iqc2
    lis_rel=block.DialogPrm(9).Data;
    if ~ischar(lis_rel)
        disp_str(9,'custom IQC relationships','char')
    end
    setvalue(pos,'iqc_custom',{'0'},{'lis_rel',lis_rel});
end

function Start(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_uncertain_cdelay';
pos=checkblockpos(block_path,glparam_name);

delay_time=block.DialogPrm(10).Data;
ext_iqc_setvalues.ext_sim_uncertain_cdelay.delay_time{pos}=delay_time;
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
    error('The number of random ~= dimension of input signal')
end

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_uncertain_cdelay';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_sim_uncertain_cdelay.sampletime{pos}(1:2)=di(1);
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
block_path=gcb;
glparam_name='ext_sim_uncertain_cdelay';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_sim_uncertain_cdelay.size{pos}(1:2)=di;
% Set compiled dimensions
block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = di;

function Outputs(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_uncertain_cdelay';
pos=checkblockpos(block_path,glparam_name);

delay_time=block.DialogPrm(10).Data;

if ~isfield(ext_iqc_setvalues.ext_sim_uncertain_cdelay,'size')
    ext_iqc_setvalues.ext_sim_uncertain_cdelay.size{pos}(1:2)=...
        block.InputPort(1).Dimensions;
end
% block.OutputPort(1).Data=zeros(block.InputPort(1).Dimensions,1);
block.OutputPort(1).Data=delay_time';
%endfunction

function Terminate(block)
iqctool('simclear')
%end Terminate

function setvalue(pos,iqc_name,out_var,other)
global ext_iqc_setvalues
str=[];
for i=1:length(out_var)
    str=[str,out_var{i},','];
end
try
    ny1=size(ext_iqc_setvalues.ext_sim_uncertain_cdelay.iqctype{pos},2);
catch err
    ny1=0;
end
if ~strcmp(str(1:end-1),'0')
    ext_iqc_setvalues.ext_sim_uncertain_cdelay.iqctype{pos}{ny1+1}(1,:)={[iqc_name,...
        '.outvariable'],...
        str};
else
    ext_iqc_setvalues.ext_sim_uncertain_cdelay.iqctype{pos}{ny1+1}(1,:)={[iqc_name,...
        '.list_iqc_rel'],...
        other{1,2}};
    other=[];
end
for i=1:size(other,1)
    
    ny2=size(ext_iqc_setvalues.ext_sim_uncertain_cdelay.iqctype{pos}{ny1+1},1);
    ext_iqc_setvalues.ext_sim_uncertain_cdelay.iqctype{pos}{ny1+1}(ny2+1,:)=...
        {[iqc_name,'.',other{i,1}],other{i,2}};
end
