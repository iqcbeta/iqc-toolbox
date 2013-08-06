function ext_iqc_popov_vect(block)
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
str='disp(''Popov'');';
set_param(gcb,'MaskDisplay',str);


function CheckParameters(block)
sign_par = block.DialogPrm(1).Data;
out_var = block.DialogPrm(2).Data;
for i=1:length(out_var)
    if ~isvarname(out_var{i})
        DAStudio.error('Simulink:block:invalidParameter');
    end
end
if length(out_var)>1,
    DAStudio.error('Simulink:block:invalidParameter');
end
switch sign_par
    case 1
        sign_par='0';
    case 2
        sign_par='+';
    case 3
        sign_par='-';
end

block_path=gcb;
ny=strfind(block_path,'/');
block_newpath=[block_path(1:ny-1),'_lin',block_path(ny:size(block_path,2))];
block_name=get_param(gcb,'name');

str=[];
for i=1:length(out_var)
    str=[str,out_var{i},','];
end

glparam_name='ext_iqc_popov_vect';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','sign_par','outvariable','path','newpath'};
alparam={block_name,sign_par,str,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);


function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_popov_vect';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_iqc_popov_vect.sampletime{pos}(1:2)=di(1);
block.InputPort(idx).Sampletime = di;
block.OutputPort(idx).Sampletime = di;

function SetOutputSampleTime(block, idx, di)


function SetInputPortDims(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_popov_vect';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_iqc_popov_vect.size{pos}(1:2)=di;
% Set compiled dimensions
block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = di;


function Outputs(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_popov_vect';
pos=checkblockpos(block_path,glparam_name);

if ~isfield(ext_iqc_setvalues.ext_iqc_popov_vect,'size')
    ext_iqc_setvalues.ext_iqc_popov_vect.size{pos}(1:2)=...
       block.InputPort(1).Dimensions;
end
%endfunction

function Terminate(block)
iqctool('simclear')
%end Terminate
