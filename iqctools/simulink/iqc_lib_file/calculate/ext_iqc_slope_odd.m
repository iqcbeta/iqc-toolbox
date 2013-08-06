function ext_iqc_slope_odd(block)
setup(block);

function setup(block)

block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 5;

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
str=['plot(0,0,100,100,[96,2],[50,50],[50,50],[95,4],',...
    '[97,80,60,40,20,4],[90,90,71,30,10,10],[81,84,60,75,60],',...
    '[72,71,71,99,71],[81,84],[69,71],[73,75],[97,99],',...
    '[75,75,75],[96,99,99]);'];
set_param(gcb,'MaskDisplay',str);


function CheckParameters(block)
pole = block.DialogPrm(1).Data;
len_def_h = block.DialogPrm(2).Data;
low_bou = block.DialogPrm(3).Data;
upp_bou = block.DialogPrm(4).Data;
out_var = block.DialogPrm(5).Data;
if ~isa(pole,'double')
    DAStudio.error('Simulink:block:invalidParameter');
end
if ~isa(len_def_h,'double') || length(len_def_h)~=1
    DAStudio.error('Simulink:block:invalidParameter');
end
if ~isa(low_bou,'double') || length(low_bou)~=1
    DAStudio.error('Simulink:block:invalidParameter');
end
if ~isa(upp_bou,'double') || length(upp_bou)~=1
    DAStudio.error('Simulink:block:invalidParameter');
end
for i=1:length(out_var)
    if ~isvarname(out_var{i})
        DAStudio.error('Simulink:block:invalidParameter');
    end
end
if length(out_var)>4,
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

glparam_name='ext_iqc_slope_odd';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','pole','len_def_h','low_bou','upp_bou',...
    'outvariable','path','newpath'};
alparam={block_name,pole,len_def_h,low_bou,upp_bou,...
    str,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_slope_odd';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_iqc_slope_odd.sampletime{pos}(1:2)=di(1);
block.InputPort(idx).Sampletime = di;
block.OutputPort(idx).Sampletime = di;

function SetOutputSampleTime(block, idx, di)


function SetInputPortDims(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_slope_odd';
pos=checkblockpos(block_path,glparam_name);

% Set compiled dimensions
block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = di;
ext_iqc_setvalues.ext_iqc_slope_odd.size{pos}(1:2)=di;

function Outputs(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_slope_odd';
pos=checkblockpos(block_path,glparam_name);

if ~isfield(ext_iqc_setvalues.ext_iqc_slope_odd,'size')
    ext_iqc_setvalues.ext_iqc_slope_odd.size{pos}(1:2)=...
        block.InputPort(1).Dimensions;
end
block.OutputPort(1).Data=zeros(block.OutputPort(1).Dimensions,1);
%endfunction

function Terminate(block)
iqctool('simclear')
%end Terminate
