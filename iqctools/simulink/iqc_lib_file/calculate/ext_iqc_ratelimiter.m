function ext_iqc_ratelimiter(block)
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
str=['plot(0,0,100,100,[18,45,45,18,18],[60,60,26,26,60],',...
    '[62,82,82,62,62],[60,60,26,26,60],[62,0],[44,44],',...
    '[79,65],[44,44],[71,71],[46,57],',...
    '[68,69,72,74,73,71,69,68,69,72,74],',...
    '[31,29,29,31,34,35,36,38,40,40,39],[31,31],[26,60],',...
    '[45,38,24],[55,55,33],[24,18],[33,33],[99,82],[44,44],',...
    '[10,12,8,10,10,90,90],[44,51,51,44,80,80,44],[7,2],[55,55]);'];
set_param(gcb,'MaskDisplay',str);

function CheckParameters(block)
pole = block.DialogPrm(1).Data;
k = block.DialogPrm(2).Data;
out_var = block.DialogPrm(3).Data;
if ~isa(pole,'double')
    DAStudio.error('Simulink:block:invalidParameter');
end
if ~isa(k,'double') || length(k)~=1
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

glparam_name='ext_iqc_ratelimiter';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','pole','k','outvariable','path','newpath'};
alparam={block_name,pole,k,str,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_ratelimiter';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_iqc_ratelimiter.sampletime{pos}(1:2)=di(1);
block.InputPort(idx).Sampletime = di;
block.OutputPort(idx).Sampletime = di;

function SetOutputSampleTime(block, idx, di)


function SetInputPortDims(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_ratelimiter';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_iqc_ratelimiter.size{pos}(1:2)=di;

block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = di;


function Outputs(block)
block_path=gcb;
glparam_name='ext_iqc_ratelimiter';
pos=checkblockpos(block_path,glparam_name);

global ext_iqc_setvalues
if ~isfield(ext_iqc_setvalues.ext_iqc_ratelimiter,'size')
    ext_iqc_setvalues.ext_iqc_ratelimiter.size{pos}(1:2)=...
       block.InputPort(1).Dimensions;
end
block.OutputPort(1).Data=zeros(block.OutputPort(1).Dimensions,1);
%endfunction

function Terminate(block)
iqctool('simclear')
%end Terminate
