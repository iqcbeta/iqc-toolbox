function ext_sim_saturation(block)
setup(block);

function setup(block)
block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 16;

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions = -1;
block.InputPort(1).DirectFeedthrough = true;
block.InputPort(1).SamplingMode = 'sample';
block.InputPort(1).SampleTime = [-1 0];

sam_time=block.DialogPrm(3).Data;
if sam_time==-2,
    sam_time=0;
end

block.OutputPort(1).Dimensions = -1;
block.OutputPort(1).SamplingMode = 'sample';
block.OutputPort(1).SampleTime = [sam_time 0];

block.RegBlockMethod('SetInputPortDimensions',@SetInputPortDims);
block.RegBlockMethod('SetInputPortSampleTime',@SetInputSampleTime);
block.RegBlockMethod('SetOutputPortSampleTime',@SetOutputSampleTime);
block.RegBlockMethod('CheckParameters',@CheckParameters);
block.RegBlockMethod('Outputs', @Outputs);
block.RegBlockMethod('Terminate', @Terminate);
block.RegBlockMethod('InitializeConditions',@InitializeConditions);

function InitializeConditions(block)
str=['x=-10:10;',...
    'for i=1:21;',...
    'if x(i)>=-3 && x(i)<=3;',...
    'y(i)=x(i);',...
    'elseif x(i)<-3;',...
    'y(i)=-3;',...
    'elseif x(i)>3;',...
    'y(i)=3;',...
    'end;',...
    'end;',...
    'x1=zeros(1,11);y1=-5:1:5;y2=zeros(1,21);plot(x,y);plot(x1,y1);plot(x,y2)'];
set_param(gcb,'MaskDisplay',str);

function CheckParameters(block)
global ext_iqc_setvalues

% Simulink Parameters Check
upp_lim=block.DialogPrm(1).Data;
if ~isa(upp_lim,'double')
    disp_str(9,'upper limit','numeric')
end
low_lim=block.DialogPrm(2).Data;
if ~isa(low_lim,'double')
    disp_str(9,'lower limit','numeric')
end
if low_lim>upp_lim
    disp_str(39,'"upper limit"','<','"lower limit"')
end

sam_tim=block.DialogPrm(3).Data;
if ~isa(sam_tim,'double')
    disp_str(9,'sample time','numeric')
end

% IQC Procedure Parameters Check
block_path=gcb;
ny=strfind(block_path,'/');
block_newpath=[block_path(1:ny-1),'_lin',block_path(ny:size(block_path,2))];
block_name=get_param(gcb,'name');

glparam_name='ext_sim_saturation';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','path','newpath'};
alparam={block_name,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);

ext_iqc_setvalues.ext_sim_saturation.iqctype{pos}='';

iqc1=block.DialogPrm(4).Data;
iqc2=block.DialogPrm(8).Data;
iqc3=block.DialogPrm(10).Data;
iqc4=block.DialogPrm(12).Data;
iqc5=block.DialogPrm(15).Data;

if ~iqc1 && ~iqc2 && ~iqc3 && ~iqc4 && ~iqc5
    DAStudio.error('Simulink:block:invalidParameter');
    iqctool('simclear')
end

if iqc1
    iqc1_pole=block.DialogPrm(6).Data;
    if ~isa(iqc1_pole,'double')
        disp_str(9,'"iqc_monotonic pole"','numeric')
    end
    iqc1_ope=block.DialogPrm(6).Data;
    iqc1_out=block.DialogPrm(7).Data;
    for i=1:length(iqc1_out)
        if ~isvarname(iqc1_out{i})
            DAStudio.error('Simulink:block:invalidParameter');
        end
    end
    switch iqc1_ope
        case 1
            ope=0;
        case 2
            ope=1;
    end
    setvalue(pos,'iqc_monotonic',iqc1_out,{'pole',iqc1_pole;...
        'ope',ope;...
        'int_dd',1});
end
if iqc2
    iqc2_out=block.DialogPrm(9).Data;
    for i=1:length(iqc2_out)
        if ~isvarname(iqc2_out{i})
            DAStudio.error('Simulink:block:invalidParameter');
        end
    end
    setvalue(pos,'iqc_sector',iqc2_out,{'low_sec',0;...
        'upp_sec',1});
end
if iqc3
    iqc3_out=block.DialogPrm(11).Data;
    for i=1:length(iqc3_out)
        if ~isvarname(iqc3_out{i})
            DAStudio.error('Simulink:block:invalidParameter');
        end
    end
    setvalue(pos,'iqc_popov',iqc3_out,{'low_sec',0;...
        'upp_sec',1;});
end
if iqc4
    iqc4_sign=block.DialogPrm(13).Data;
    iqc4_out=block.DialogPrm(14).Data;
    for i=1:length(iqc4_out)
        if ~isvarname(iqc4_out{i})
            DAStudio.error('Simulink:block:invalidParameter');
        end
    end
    switch iqc4_sign
        case 1
            iqc4_sign='0';
        case 2
            iqc4_sign='+';
        case 3
            iqc4_sign='-';
    end
    setvalue(pos,'iqc_popov_vect',iqc4_out,{'sign_par',iqc4_sign});
end
if iqc5
    lis_rel=block.DialogPrm(16).Data;
    if ~ischar(lis_rel)
        disp_str(9,'custom IQC relationships','char')
    end
    setvalue(pos,'iqc_custom',{'0'},{'lis_rel',lis_rel});
end

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_saturation';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_sim_saturation.sampletime{pos}(1:2)=di(1);
if di(1)==0,
    block.InputPort(idx).Sampletime = [0 di(2)];
    block.OutputPort(idx).Sampletime = [0 0];
else
    block.InputPort(idx).Sampletime = di;
    block.OutputPort(idx).Sampletime = [di(1) 0];
end

function SetOutputSampleTime(block, idx, di)


function SetInputPortDims(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_saturation';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_sim_saturation.size{pos}(1:2)=di;
% Set compiled dimensions
block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = di;


function Outputs(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_saturation';
pos=checkblockpos(block_path,glparam_name);

if ~isfield(ext_iqc_setvalues.ext_sim_saturation,'size')
    ext_iqc_setvalues.ext_sim_saturation.size{pos}(1:2)=...
        block.InputPort(1).Dimensions;
end

upp_lim=block.DialogPrm(1).Data;
low_lim=block.DialogPrm(2).Data;

ny=size(block.InputPort(1).Data,1);

for i=1:ny,
    
    output_value=block.InputPort(1).Data(i,1);
    if output_value<upp_lim && low_lim<output_value
        block.OutputPort(1).Data(i,1)=output_value;
    elseif upp_lim<=output_value
        block.OutputPort(1).Data(i,1)=upp_lim;
    elseif output_value<=low_lim
        block.OutputPort(1).Data(i,1)=low_lim;
    end
end

function Terminate(block)
iqctool('simclear')

function setvalue(i_pos,iqc_name,out_var,other)
global ext_iqc_setvalues
str=[];
for i=1:length(out_var)
    str=[str,out_var{i},','];
end
try
    ny1=size(ext_iqc_setvalues.ext_sim_saturation.iqctype{i_pos},2);
catch err
    ny1=0;
end
if ~strcmp(str(1:end-1),'0')
    ext_iqc_setvalues.ext_sim_saturation.iqctype{i_pos}{ny1+1}(1,:)={[iqc_name,...
        '.outvariable'],...
        str};
else
    ext_iqc_setvalues.ext_sim_saturation.iqctype{i_pos}{ny1+1}(1,:)={[iqc_name,...
        '.list_iqc_rel'],...
        other{1,2}};
    other=[];
end
for i=1:size(other,1)
    
    ny2=size(ext_iqc_setvalues.ext_sim_saturation.iqctype{i_pos}{ny1+1},1);
    ext_iqc_setvalues.ext_sim_saturation.iqctype{i_pos}{ny1+1}(ny2+1,:)=...
        {[iqc_name,'.',other{i,1}],other{i,2}};
end
