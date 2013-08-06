function ext_iqc_dzn_e_odd(block)
setup(block);

function setup(block)

block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 4;

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
str=['plot(0,0,100,100,[60,30,30,60,60],[90,90,60,60,90],',...
    '[30,30,60,60,30],[40,10,10,40,40],[30,10,10,3,30],',...
    '[25,25,76,76,76],[80,80,60],[72,25,25],[60,77],[76,76],',...
    '[82,80,77,80,82,98],[76,72,76,80,76,76],[58,50,40,32],',...
    '[88,76,76,63],[42,45,45,50],[15,17,34,36]);'];
set_param(gcb,'MaskDisplay',str);


function CheckParameters(block)
pole = block.DialogPrm(1).Data;
def_num_h = block.DialogPrm(2).Data;
slo_bou = block.DialogPrm(3).Data;
out_var = block.DialogPrm(4).Data;
if ~isa(pole,'double')
    DAStudio.error('Simulink:block:invalidParameter');
end
if ~isa(def_num_h,'double')|| length(def_num_h)~=1
    DAStudio.error('Simulink:block:invalidParameter');
end
if ~isa(slo_bou,'double') || length(slo_bou)~=1
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

glparam_name='ext_iqc_dzn_e_odd';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','pole','def_num_h','slo_bou',...
    'outvariable','path','newpath'};
alparam={block_name,pole,def_num_h,slo_bou,...
    str,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_dzn_e_odd';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_iqc_dzn_e_odd.sampletime{pos}(1:2)=di(1);
block.InputPort(idx).Sampletime = di;
block.OutputPort(idx).Sampletime = di;

function SetOutputSampleTime(block, idx, di)


function SetInputPortDims(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_dzn_e_odd';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_iqc_dzn_e_odd.size{pos}(1:2)=di;
% Set compiled dimensions
block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = di;


function Outputs(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_dzn_e_odd';
pos=checkblockpos(block_path,glparam_name);
if ~isfield(ext_iqc_setvalues.ext_iqc_dzn_e_odd,'size')
    ext_iqc_setvalues.ext_iqc_dzn_e_odd.size{pos}(1:2)=...
        block.InputPort(1).Dimensions;
end
%endfunction

function Terminate(block)
iqctool('simclear')
%end Terminate
