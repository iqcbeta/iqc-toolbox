function ext_iqc_polytope_stvp(block)
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
str='disp(''STV\nPolytopic'');';
set_param(gcb,'MaskDisplay',str);


function CheckParameters(block)
arr_del = block.DialogPrm(1).Data;
arr_ddel = block.DialogPrm(2).Data;
str_mat = block.DialogPrm(3).Data;
out_var = block.DialogPrm(4).Data;
for i=1:length(out_var)
    if ~isvarname(out_var{i})
        DAStudio.error('Simulink:block:invalidParameter');
    end
end
if isempty(arr_del),
    DAStudio.error('Simulink:block:invalidParameter');
end
if isempty(arr_ddel),
    DAStudio.error('Simulink:block:invalidParameter');
end
% if isempty(str_mat),
%     DAStudio.error('Simulink:block:invalidParameter');
% end
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

glparam_name='ext_iqc_polytope_stvp';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','arr_del','arr_ddel','str_mat','outvariable','path','newpath'};
alparam={block_name,arr_del,arr_ddel,str_mat,str,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_polytope_stvp';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_iqc_polytope_stvp.sampletime{pos}(1:2)=di(1);
block.InputPort(idx).Sampletime = di;
block.OutputPort(idx).Sampletime = di;

function SetOutputSampleTime(block, idx, di)


function SetInputPortDims(block, idx, di)
arr_del = block.DialogPrm(1).Data;
ny=size(arr_del{1},1);
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_polytope_stvp';
pos=checkblockpos(block_path,glparam_name);

% Set compiled dimensions
block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = ny;
ext_iqc_setvalues.ext_iqc_polytope_stvp.size{pos}(1)=di;
ext_iqc_setvalues.ext_iqc_polytope_stvp.size{pos}(2)=ny;

function Outputs(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_iqc_polytope_stvp';
pos=checkblockpos(block_path,glparam_name);

if ~isfield(ext_iqc_setvalues.ext_iqc_polytope_stvp,'size')
    out_dim = block.DialogPrm(1).Data;
    ext_iqc_setvalues.ext_iqc_polytope_stvp.size{pos}(1)=...
        block.InputPort(1).Dimensions;
    if out_dim==-1,
        block.OutputPort(1).Dimensions =  block.InputPort(1).Dimensions;
        ext_iqc_setvalues.ext_iqc_polytope_stvp.size{pos}(2)=...
            block.InputPort(1).Dimensions;
    else
        block.OutputPort(1).Dimensions = out_dim;
        ext_iqc_setvalues.ext_iqc_polytope_stvp.size{pos}(2)=out_dim;
    end
end
block.OutputPort(1).Data=zeros(block.OutputPort(1).Dimensions,1);
%endfunction

function Terminate(block)
iqctool('simclear')
%end Terminate
