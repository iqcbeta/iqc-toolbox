function ext_iqc_popov(block)
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
str=['plot(0,0,100,100,[100,-0],[54,54],[50,50],[0,100],[66,32],',...
    '[100,0],[-0,100],[37,71],[20,-0],[0,32],[33,13],[5,41],',...
    '[40,28],[24,46],[47,42],[44,51],[65,57],[60,73],[78,63],',...
    '[64,90],[95,77],[69,100],[91,85,82,78,77,75,74,66,61,56,',...
    '54,50,48,45,43,42,40,38,35,31,27,24,21,17,14,11,10,8,7],',...
    '[96,96,93,90,82,76,71,66,66,63,60,54,51,49,46,42,38,34,32,',...
    '32,32,32,32,31,29,27,24,18,13]);'];
set_param(gcb,'MaskDisplay',str);

function CheckParameters(block)
low_sec = block.DialogPrm(1).Data;
upp_sec = block.DialogPrm(2).Data;
out_var = block.DialogPrm(3).Data;
if ~isa(low_sec,'double') || length(low_sec)~=1
    DAStudio.error('Simulink:block:invalidParameter');
end
if ~isa(upp_sec,'double') || length(upp_sec)~=1
    DAStudio.error('Simulink:block:invalidParameter');
end
if upp_sec<low_sec
    DAStudio.error('Simulink:block:invalidParameter');
end
for i=1:length(out_var)
    if ~isvarname(out_var{i})
        DAStudio.error('Simulink:block:invalidParameter');
    end
end
if length(out_var)>2,
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

glparam_name='ext_iqc_popov';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','low_sec','upp_sec','outvariable','path','newpath'};
alparam={block_name,low_sec,upp_sec,str,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_popov';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_iqc_popov.sampletime{pos}(1:2)=di(1);
block.InputPort(idx).Sampletime = di;
block.OutputPort(idx).Sampletime = di;

function SetOutputSampleTime(block, idx, di)


function SetInputPortDims(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_popov';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_iqc_popov.size{pos}(1:2)=di;
% Set compiled dimensions
block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = di;


function Outputs(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_popov';
pos=checkblockpos(block_path,glparam_name);

if ~isfield(ext_iqc_setvalues.ext_iqc_popov,'size')
    ext_iqc_setvalues.ext_iqc_popov.size{pos}(1:2)=...
        block.InputPort(1).Dimensions;
end
block.OutputPort(1).Data=zeros(block.OutputPort(1).Dimensions,1);
%endfunction

function Terminate(block)
iqctool('simclear')
%end Terminate
