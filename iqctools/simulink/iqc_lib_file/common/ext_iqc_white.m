function ext_iqc_white(block)
setup(block);

function setup(block)

block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 5;

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions = -1;
block.InputPort(1).DirectFeedthrough = true;
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
block.RegBlockMethod('Terminate', @Terminate);

function CheckParameters(block)
sig_size = block.DialogPrm(1).Data;
ban_whi = block.DialogPrm(2).Data;
pole = block.DialogPrm(3).Data;
out_var = block.DialogPrm(4).Data;
for i=1:length(out_var)
    if ~isvarname(out_var{i})
        DAStudio.error('Simulink:block:invalidParameter');
    end
end
if ~isa(sig_size,'double') || length(sig_size)~=1 || sig_size<=0
    DAStudio.error('Simulink:block:invalidParameter');
end
if ~isa(ban_whi,'double') || length(ban_whi)~=1 || ban_whi<=0
    DAStudio.error('Simulink:block:invalidParameter');
end
if ~isa(pole,'double')
    DAStudio.error('Simulink:block:invalidParameter');
end
if length(out_var)>2,
    DAStudio.error('Simulink:block:invalidParameter');
end

block_path=gcb;
dot_pos=strfind(block_path,'/');
block_path=block_path(1:dot_pos(length(dot_pos))-1);
glparam_name='ext_iqc_white';
pos=checkblockpos(block_path,glparam_name);

block_name=get_param(block_path,'name');
dot_pos=strfind(block_path,'/');
block_newpath=[block_path(1:dot_pos(1)-1),'_lin',...
    block_path(dot_pos(1):length(block_path))];

str=[];
for i=1:length(out_var)
    str=[str,out_var{i},','];
end

alparam_name={'name','ban_whi','pole',...
    'outvariable','path','newpath'};
alparam={block_name,ban_whi,pole,...
    str,block_path,block_newpath};

setblockparam(pos,glparam_name,alparam_name,alparam);

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
sam_time = block.DialogPrm(5).Data;
block_path=gcb;
dot_pos=strfind(block_path,'/');
block_path=block_path(1:dot_pos(length(dot_pos))-1);
glparam_name='ext_iqc_white';
pos=checkblockpos(block_path,glparam_name);
ext_iqc_setvalues.ext_iqc_white.sampletime{pos}(1:2)=sam_time;
block.InputPort(idx).Sampletime = di;
block.OutputPort(idx).Sampletime = [sam_time 0];

function SetOutputSampleTime(block, idx, di)

function SetInputPortDims(block, idx, di)
global ext_iqc_setvalues
sig_size = block.DialogPrm(1).Data;
block_path=gcb;
dot_pos=strfind(block_path,'/');
block_path=block_path(1:dot_pos(length(dot_pos))-1);
glparam_name='ext_iqc_white';
pos=checkblockpos(block_path,glparam_name);
ext_iqc_setvalues.ext_iqc_white.sig_size{pos}(1)=0;
ext_iqc_setvalues.ext_iqc_white.sig_size{pos}(2)=sig_size;
block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = sig_size;

function Outputs(block)
sig_size = block.DialogPrm(1).Data;
block.OutputPort(1).Data=ones(sig_size,1)*block.InputPort(1).Data;

function Terminate(block)
iqctool('simclear')
