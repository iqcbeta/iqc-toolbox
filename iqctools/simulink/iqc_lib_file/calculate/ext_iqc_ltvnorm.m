function ext_iqc_ltvnorm(block)
setup(block);

function setup(block)

block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 3;

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
block.RegBlockMethod('Terminate', @Terminate);
block.RegBlockMethod('InitializeConditions',@InitializeConditions);

function InitializeConditions(block)
str='disp(''\fontsize{16}||\Delta|| < k'', ''texmode'',''on'');';
set_param(gcb,'MaskDisplay',str);

function CheckParameters(block)
out_dim = block.DialogPrm(1).Data;
max_amp = block.DialogPrm(2).Data;
out_var = block.DialogPrm(3).Data;
if ~isa(out_dim,'double') || length(out_dim)~=1 || out_dim<1
    DAStudio.error('Simulink:block:invalidParameter');
end
if ~isa(max_amp,'double') || length(max_amp)~=1 || max_amp<=0
    DAStudio.error('Simulink:block:invalidParameter');
end
for i=1:length(out_var)
    if ~isvarname(out_var{i})
        DAStudio.error('Simulink:block:invalidParameter');
    end
end

if length(out_var)>1,
    DAStudio.error('Simulink:block:invalidParameter');
end

block_path=gcb;
ny=strfind(block_path,'/');
block_newpath=[block_path(1:ny-1),'_lin',block_path(ny:size(block_path,2))];
block_name=get_param(gcb,'name');

str=[];
for i=1:length(out_var)
    str=[str,out_var{i},','];
end

glparam_name='ext_iqc_ltvnorm';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','max_amp','outvariable','path','newpath'};
alparam={block_name,max_amp,str,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_ltvnorm';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_iqc_ltvnorm.sampletime{pos}(1:2)=di(1);
block.InputPort(idx).Sampletime = di;
block.OutputPort(idx).Sampletime = di;

function SetOutputSampleTime(block, idx, di)

function SetInputPortDims(block, idx, di)
global ext_iqc_setvalues
out_dim=block.DialogPrm(1).Data;

block_path=gcb;
glparam_name='ext_iqc_ltvnorm';
pos=checkblockpos(block_path,glparam_name);

block.InputPort(idx).Dimensions = di;
ext_iqc_setvalues.ext_iqc_ltvnorm.size{pos}(1)=di;
if out_dim==-1,
    block.OutputPort(idx).Dimensions = di;
    ext_iqc_setvalues.ext_iqc_ltvnorm.out_dim{pos}=di;
    ext_iqc_setvalues.ext_iqc_ltvnorm.size{pos}(2)=di;
else
    block.OutputPort(idx).Dimensions = out_dim;
    ext_iqc_setvalues.ext_iqc_ltvnorm.out_dim{pos}=out_dim;
    ext_iqc_setvalues.ext_iqc_ltvnorm.size{pos}(2)=out_dim;
end

function Outputs(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_ltvnorm';
pos=checkblockpos(block_path,glparam_name);

if ~isfield(ext_iqc_setvalues.ext_iqc_ltvnorm,'size')
    out_dim = block.DialogPrm(1).Data;
    ext_iqc_setvalues.ext_iqc_ltvnorm.size{pos}(1)=...
        block.InputPort(1).Dimensions;
    if out_dim==-1,
        block.OutputPort(1).Dimensions =  block.InputPort(1).Dimensions;
        ext_iqc_setvalues.ext_iqc_ltvnorm.size{pos}(2)=...
            block.InputPort(1).Dimensions;
    else
        block.OutputPort(1).Dimensions = out_dim;
        ext_iqc_setvalues.ext_iqc_ltvnorm.size{pos}(2)=out_dim;
    end
end
block.OutputPort(1).Data=zeros(block.OutputPort(1).Dimensions,1);

function Terminate(block)
iqctool('simclear')
