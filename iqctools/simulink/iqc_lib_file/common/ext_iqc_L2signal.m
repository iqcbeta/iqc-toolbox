function ext_iqc_L2signal(block)
setup(block);

function setup(block)

sig_size=block.DialogPrm(1).Data;
sam_time=block.DialogPrm(2).Data;
sig_type=block.DialogPrm(3).Data;
sig_fre=block.DialogPrm(6).Data;

if sig_type==6,
    block.NumInputPorts  = 1;
    block.InputPort(1).DimensionsMode = 'Fixed';
    block.InputPort(1).Dimensions  = -1;
    block.InputPort(1).SamplingMode = 'sample';
else
    block.NumInputPorts  = 0;
end

block.NumOutputPorts = 1;
block.SetPreCompPortInfoToDefaults;
block.OutputPort(1).DimensionsMode = 'Fixed';
block.OutputPort(1).Dimensions  = sig_size;
block.OutputPort(1).SamplingMode = 'sample';

block.NumDialogPrms = 7;

switch sig_type
    case 1
        if sam_time==-2 || sam_time==0,
            sam_time=0.1;
        end
    case {2,3}
        if sam_time==-2 || sam_time==0,
            sam_time=0;
        end
    case {4,5,6}
        if sam_time==-2 || sam_time==0,
            sam_time=0;
        end
end

block.OutputPort(1).SampleTime = [sam_time 0];
block.SimStateCompliance = 'DefaultSimState';

block.RegBlockMethod('Outputs', @Outputs);
block.RegBlockMethod('Terminate', @Terminate);
block.RegBlockMethod('SetInputPortDimensions',@SetInputPortDims);
block.RegBlockMethod('SetInputPortSampleTime',@SetInputSampleTime);
block.RegBlockMethod('SetOutputPortDimensions',@SetOutputPortDims);
block.RegBlockMethod('SetOutputPortSampleTime',@SetOutputSampleTime);
block.RegBlockMethod('Start',@Start);
block.RegBlockMethod('InitializeConditions',@InitializeConditions);

function InitializeConditions(block)
str='plot(0:0.1:10,[sin(0:0.1:6.2) zeros(1,38)])';
set_param(gcb,'MaskDisplay',str);

function Start(block)
global ext_iqc_setvalues
sig_size=block.DialogPrm(1).Data;
sam_time=block.DialogPrm(2).Data;
block_path=gcb;
ny=strfind(block_path,'/');
block_newpath=[block_path(1:ny-1),'_lin',block_path(ny:size(block_path,2))];
block_name=get_param(gcb,'name');
glparam_name='ext_iqc_L2signal';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','path','newpath'};
alparam={block_name,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);
ext_iqc_setvalues.ext_iqc_L2signal.size{pos}(1:2)=[0 sig_size];
ext_iqc_setvalues.ext_iqc_L2signal.sampletime{pos}(1:2)=[sam_time sam_time];

function SetOutputSampleTime(block, idx, di)

function SetOutputPortDims(block, idx, di)

function SetInputSampleTime(block, idx, di)
block.InputPort(idx).SampleTime=di;

function SetInputPortDims(block, idx, di)
block.InputPort(idx).Dimensions  = di;

function Outputs(block)
sig_size=block.DialogPrm(1).Data;
sig_type=block.DialogPrm(3).Data;
tru_time=block.DialogPrm(4).Data;
sig_amp=block.DialogPrm(5).Data;
sig_fre=block.DialogPrm(6).Data;
user_cus=block.DialogPrm(7).Data;
t=block.CurrentTime;
if t<tru_time,
    switch sig_type,
        case 1
            block.OutputPort(1).Data=sig_amp*ones(sig_size,1);
        case 2
            test_time=fix(sig_fre*t);
            if ~mod(test_time,2)
                block.OutputPort(1).Data=sig_amp*ones(sig_size,1);
            else
                block.OutputPort(1).Data=zeros(sig_size,1);
            end
        case 3
            block.OutputPort(1).Data=ones(sig_size,1)*sig_amp*sin(sig_fre*t);
        case 4
            block.OutputPort(1).Data=ones(sig_size,1)*eval(user_cus);
        case 5
            block.OutputPort(1).Data=ones(sig_size,1)*sig_amp*rand(1);
        case 6
            block.OutputPort(1).Data=ones(sig_size,1)*block.InputPort(1).Data;
    end
else
    block.OutputPort(1).Data=zeros(sig_size,1);
end

function Terminate(block)
iqctool('simclear')
