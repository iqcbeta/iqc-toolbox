function ext_iqc_estimation(block)

setup(block);

%endfunction

function setup(block)

block.NumInputPorts  = 5;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 5;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions = -1;
block.InputPort(1).DirectFeedthrough = false;
block.InputPort(1).SamplingMode = 'sample';
block.InputPort(1).SampleTime = [-1 0];
block.InputPort(2).Dimensions = -1;
block.InputPort(2).DirectFeedthrough = false;
block.InputPort(2).SampleTime = [-1 0];
block.InputPort(2).SamplingMode = 'sample';
block.InputPort(3).Dimensions = -1;
block.InputPort(3).DirectFeedthrough = false;
block.InputPort(3).SampleTime = [-1 0];
block.InputPort(3).SamplingMode = 'sample';
block.InputPort(4).Dimensions = -1;
block.InputPort(4).DirectFeedthrough = false;
block.InputPort(4).SampleTime = [-1 0];
block.InputPort(4).SamplingMode = 'sample';
block.InputPort(5).Dimensions = -1;
block.InputPort(5).DirectFeedthrough = false;
block.InputPort(5).SampleTime = [-1 0];
block.InputPort(5).SamplingMode = 'sample';
block.OutputPort(1).Dimensions = -1;
block.OutputPort(1).SampleTime = [-1 0];
block.OutputPort(1).SamplingMode = 'sample';

% block.SimStateCompliance = 'DefaultSimState';
block.RegBlockMethod('SetInputPortDimensions',@SetInputPortDims);
block.RegBlockMethod('SetInputPortSampleTime',@SetInputSampleTime);
block.RegBlockMethod('SetOutputPortSampleTime',@SetOutputSampleTime);
block.RegBlockMethod('Outputs', @Outputs);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Terminate', @Terminate);
% block.RegBlockMethod('InitializeConditions',@InitializeConditions);
block.RegBlockMethod('CheckParameters',@CheckParameters);
block.RegBlockMethod('Update', @Update);
block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('PostPropagationSetup', @PostPropagationSetup);

function PostPropagationSetup(block)
[gain,design_sys]=iqc_gui(gcs);

desmet=block.DialogPrm(1).Data;
switch desmet
    case 1
        
end


function CheckParameters(block)
desmet=block.DialogPrm(1).Data;
save_lin_com=block.DialogPrm(2).Data;
uselmitool=block.DialogPrm(3).Data;
solver=block.DialogPrm(4).Data;
lmi_par=block.DialogPrm(5).Data;
i_pos=set_init(desmet,save_lin_com,uselmitool,solver,lmi_par);

% function Start(block)
% global ext_iqc_setvalues
% per_type=block.DialogPrm(2).Data;
% save_lin_com=block.DialogPrm(3).Data;
% uselmitool=block.DialogPrm(6).Data;
% solver=block.DialogPrm(7).Data;
% lmi_par=block.DialogPrm(8).Data;
% i_pos=set_init(per_type,save_lin_com,uselmitool,solver,lmi_par);
% ext_iqc_setvalues.ext_iqc_estimation.size{i_pos}(1)=...
%     block.InputPort(2).Dimensions;
% ext_iqc_setvalues.ext_iqc_estimation.sampletime{i_pos}(1)=...
%     block.InputPort(2).Sampletime(1);
% ext_iqc_setvalues.ext_iqc_estimation.size{i_pos}(2)=...
%     block.InputPort(1).Dimensions;
% ext_iqc_setvalues.ext_iqc_estimation.sampletime{i_pos}(2)=...
%     block.InputPort(1).Sampletime(1);

function SetInputSampleTime(block, idx, di)
if idx==1
    global ext_iqc_setvalues %#ok<TLEV>
    i_pos=find_ipos(gcb);
    ext_iqc_setvalues.ext_iqc_estimation.sampletime{i_pos}(1:2)=di(1);
    block.OutputPort(idx).Sampletime = di;
end
block.InputPort(idx).Sampletime = di;

function SetOutputSampleTime(block, idx, di)

function SetInputPortDims(block, idx, di)
global ext_iqc_setvalues
block.InputPort(idx).Dimensions = di;
ext_iqc_setvalues.ext_iqc_estimation.size{i_pos}(idx)=di;
if idx==3
block.OutputPort(idx).Dimensions = di;
end

function Outputs(block)
% if ext_iqc_setvalues.ext_iqc_estimation.isrun~=1
%     ext_iqc_setvalues.ext_iqc_estimation.isrun=1;
% end
%endfunction

function Terminate(block)



function i_pos=set_init(desmet,save_lin_com,uselmitool,solver,lmi_par)
global ext_iqc_setvalues
block_path=gcb;
try
    if isempty(ext_iqc_setvalues.ext_iqc_estimation.idx)
        i_pos=0;
    else
        i_pos=ext_iqc_setvalues.ext_iqc_estimation.idx;
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

if i_pos>1
    error('The IQC Estimation Block no more than one. '),
    clear global ext_iqc_setvalues
end
% ext_iqc_setvalues.ext_iqc_estimation.isrun{i_pos}=0;
ext_iqc_setvalues.ext_iqc_estimation.save_lin_com{i_pos}=save_lin_com;
ext_iqc_setvalues.ext_iqc_estimation.path{i_pos}=block_path;
ny=strfind(block_path,'/');
block_path=[block_path(1:ny-1),'_lin',block_path(ny:end)];
ext_iqc_setvalues.ext_iqc_estimation.newpath{i_pos}=block_path;
ext_iqc_setvalues.ext_iqc_estimation.name{i_pos}=block_name;
switch desmet
    case 1
        method ='l2ctCW';
    case 2
        method='l2ctKC';
    case 3
        method='l2ctKC2';
    case 4
        method='l2dtKC';
    case 5
        method='l2dtKC2';
end
ext_iqc_setvalues.ext_iqc_estimation.method{i_pos}=method;
ext_iqc_setvalues.ext_iqc_estimation.uselmitool{i_pos}=uselmitool;
ext_iqc_setvalues.ext_iqc_estimation.solver{i_pos}=solver;
if ~iscell(lmi_par)
    lmi_par={lmi_par};
end
ext_iqc_setvalues.ext_iqc_estimation.lmi_par{i_pos}=lmi_par;
ext_iqc_setvalues.ext_iqc_estimation.idx=i_pos;
ext_iqc_setvalues.ext_iqc_estimation.jump_idx=0;

function findi=find_ipos(gcb)
global ext_iqc_setvalues
i_pos=ext_iqc_setvalues.ext_iqc_estimation.idx;
for i=1:i_pos,
    if strcmp(ext_iqc_setvalues.ext_iqc_estimation.path(i),gcb)
        findi=i;
        return
    end
end
findi=[];