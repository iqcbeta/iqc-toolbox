function ext_iqc_cdelay(block)

setup(block);

%endfunction

function setup(block)

block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 3;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% block.InputPort(1).DimensionsMode = 'Fixed';
block.InputPort(1).Dimensions = -1;
block.InputPort(1).DirectFeedthrough = false;
block.InputPort(1).SamplingMode = 'sample';
block.InputPort(1).SampleTime = [-1 0];

% block.OutputPort(1).DimensionsMode = 'Fixed';
block.OutputPort(1).Dimensions = -1;
block.OutputPort(1).SamplingMode = 'sample';
block.OutputPort(1).SampleTime = [-1 0];

% % Register sample times
% block.SampleTimes = [-1 0];

% block.SimStateCompliance = 'DefaultSimState';
% block.RegBlockMethod('InitializeCondition', @InitializeCondition);
block.RegBlockMethod('SetInputPortDimensions',@SetInputPortDims);
block.RegBlockMethod('SetInputPortSampleTime',@SetInputSampleTime);
block.RegBlockMethod('SetOutputPortSampleTime',@SetOutputSampleTime);
block.RegBlockMethod('CheckParameters',@CheckParameters);
block.RegBlockMethod('Outputs', @Outputs);
block.RegBlockMethod('Terminate', @Terminate);


% function InitializeCondition(block)

function CheckParameters(block)
max_tim = block.DialogPrm(1).Data;
pole = block.DialogPrm(2).Data;
out_var = block.DialogPrm(3).Data;
if ~strcmp(class(max_tim), 'double')
    DAStudio.error('Simulink:block:invalidParameter');
end
if ~strcmp(class(pole), 'double')
    DAStudio.error('Simulink:block:invalidParameter');
end
for i=1:length(out_var)
    if ~isvarname(out_var{i})
        DAStudio.error('Simulink:block:invalidParameter');
    end
end
if length(max_tim)~=1,
    DAStudio.error('Simulink:block:invalidParameter');
end
if length(out_var)>2,
    DAStudio.error('Simulink:block:invalidParameter');
end
set_init(max_tim,pole,out_var);

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
i_pos=find_ipos(gcb);
ext_iqc_setvalues.ext_iqc_cdelay.sampletime{i_pos}(1:2)=di(1);
block.InputPort(idx).Sampletime = di;
block.OutputPort(idx).Sampletime = di;

function SetOutputSampleTime(block, idx, di)


function SetInputPortDims(block, idx, di)
global ext_iqc_setvalues
i_pos=find_ipos(gcb);
ext_iqc_setvalues.ext_iqc_cdelay.size{i_pos}(1:2)=di;
% Set compiled dimensions
block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = di;


function Outputs(block)
i_pos=find_ipos(gcb);
global ext_iqc_setvalues
if ~isfield(ext_iqc_setvalues.ext_iqc_cdelay,'size')
    ext_iqc_setvalues.ext_iqc_cdelay.size{i_pos}(1:2)=...
       block.InputPort(1).Dimensions;
end
block.OutputPort(1).Data=zeros(block.OutputPort(1).Dimensions,1);
%endfunction

function Terminate(block)
iqctool('clear')
%end Terminate

function set_init(max_tim,pole,out_var)
global ext_iqc_setvalues
block_path=gcb;
try
    if isempty(ext_iqc_setvalues.ext_iqc_cdelay.idx)
        i_pos=0;
    else
        i_pos=ext_iqc_setvalues.ext_iqc_cdelay.idx;
        try
            findi=find_ipos(gcb);
            if ~isempty(findi),
                return
            end
        catch err
        end
    end
catch err
    i_pos=0;
end
block_name=get_param(gcb,'name');
i_pos=i_pos+1;
str=[];
for i=1:length(out_var)
    str=[str,out_var{i},','];
end

ext_iqc_setvalues.ext_iqc_cdelay.path{i_pos}=block_path;
ny=strfind(block_path,'/');
block_path=[block_path(1:ny-1),'_lin',block_path(ny:end)];
ext_iqc_setvalues.ext_iqc_cdelay.newpath{i_pos}=block_path;
ext_iqc_setvalues.ext_iqc_cdelay.name{i_pos}=block_name;
ext_iqc_setvalues.ext_iqc_cdelay.idx=i_pos;
ext_iqc_setvalues.ext_iqc_cdelay.max_tim{i_pos}=max_tim;
ext_iqc_setvalues.ext_iqc_cdelay.pole{i_pos}=pole;
ext_iqc_setvalues.ext_iqc_cdelay.outvariable{i_pos}=str;

function findi=find_ipos(gcb)
global ext_iqc_setvalues
i_pos=ext_iqc_setvalues.ext_iqc_cdelay.idx;
for i=1:i_pos,
    if strcmp(ext_iqc_setvalues.ext_iqc_cdelay.path(i),gcb)
        findi=i;
        return
    end
end
findi=[];