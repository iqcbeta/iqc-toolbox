function ext_iqc_custom(block)
setup(block);

function setup(block)

block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 2;

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
str='disp(''Custom'')';
set_param(gcb,'MaskDisplay',str);

function CheckParameters(block)
out_dim = block.DialogPrm(1).Data;
list_iqc_rel = block.DialogPrm(2).Data;
if ~isa(list_iqc_rel,'char')
    DAStudio.error('Simulink:block:invalidParameter');
end
if ~isa(out_dim,'double') || length(out_dim)~=1
    DAStudio.error('Simulink:block:invalidParameter');
end

block_path=gcb;
ny=strfind(block_path,'/');
block_newpath=[block_path(1:ny-1),'_lin',block_path(ny:size(block_path,2))];
block_name=get_param(gcb,'name');

glparam_name='ext_iqc_custom';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','list_iqc_rel','path','newpath'};
alparam={block_name,list_iqc_rel,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);


function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_custom';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_iqc_custom.sampletime{pos}(1:2)=di(1);
block.InputPort(idx).Sampletime = di;
block.OutputPort(idx).Sampletime = di;

function SetOutputSampleTime(block, idx, di)


function SetInputPortDims(block, idx, di)
out_dim = block.DialogPrm(1).Data;
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_custom';
pos=checkblockpos(block_path,glparam_name);

if out_dim==-1,
    ext_iqc_setvalues.ext_iqc_custom.size{pos}(1:2)=di;
    % Set compiled dimensions
    block.InputPort(idx).Dimensions = di;
    block.OutputPort(idx).Dimensions = di;
else
    ext_iqc_setvalues.ext_iqc_custom.size{pos}(1)=di;
    ext_iqc_setvalues.ext_iqc_custom.size{pos}(2)=out_dim;
    % Set compiled dimensions
    block.InputPort(idx).Dimensions = di;
    block.OutputPort(idx).Dimensions = out_dim;
end

function Outputs(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_custom';
pos=checkblockpos(block_path,glparam_name);

if ~isfield(ext_iqc_setvalues.ext_iqc_custom,'size')
    out_dim = block.DialogPrm(1).Data;
    ext_iqc_setvalues.ext_iqc_custom.size{pos}(1)=...
        block.InputPort(1).Dimensions;
    if out_dim==-1,
        block.OutputPort(1).Dimensions =  block.InputPort(1).Dimensions;
        ext_iqc_setvalues.ext_iqc_custom.size{pos}(2)=...
            block.InputPort(1).Dimensions;
    else
        block.OutputPort(1).Dimensions = out_dim;
        ext_iqc_setvalues.ext_iqc_custom.size{pos}(2)=out_dim;
    end
end
block.OutputPort(1).Data=zeros(block.OutputPort(1).Dimensions,1);
%endfunction

function Terminate(block)
iqctool('simclear')
%end Terminate
