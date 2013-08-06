function ext_iqc_ltigain(block)
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
str='disp(''\fontsize{16}|d| < 1'', ''texmode'',''on'');';
set_param(gcb,'MaskDisplay',str);

function CheckParameters(block)
pole = block.DialogPrm(1).Data;
out_var = block.DialogPrm(2).Data;
if ~isa(pole,'double')
    DAStudio.error('Simulink:block:invalidParameter');
end
for i=1:length(out_var)
    if ~isvarname(out_var{i})
        DAStudio.error('Simulink:block:invalidParameter');
    end
end
if length(out_var)>5,
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

glparam_name='ext_iqc_ltigain';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','pole','outvariable','path','newpath'};
alparam={block_name,pole,str,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_ltigain';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_iqc_ltigain.sampletime{pos}(1:2)=di(1);
block.InputPort(idx).Sampletime = di;
block.OutputPort(idx).Sampletime = di;
sig_dim=block.InputPort(1).Dimensions;
ext_iqc_setvalues.ext_iqc_ltigain.size{pos}(1:2)=[sig_dim sig_dim];

function SetOutputSampleTime(block, idx, di)

function SetOutputPortDims(block, idx, di)

function SetInputPortDims(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_ltigain';
pos=checkblockpos(block_path,glparam_name);
ext_iqc_setvalues.ext_iqc_ltigain.size{pos}(1:2)=di;

% Set compiled dimensions
block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = di;

function Outputs(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_ltigain';
pos=checkblockpos(block_path,glparam_name);

if ~isfield(ext_iqc_setvalues.ext_iqc_ltigain,'size')
    ext_iqc_setvalues.ext_iqc_ltigain.size{pos}(1:2)=...
       block.InputPort(1).Dimensions;
end
block.OutputPort(1).Data=zeros(block.OutputPort(1).Dimensions,1);
%endfunction

function Terminate(block)
iqctool('simclear')
%end Terminate
