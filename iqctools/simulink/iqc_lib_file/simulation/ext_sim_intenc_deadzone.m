function ext_sim_intenc_deadzone(block)
setup(block);

function setup(block)

block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 13;

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions = -1;
block.InputPort(1).DirectFeedthrough = true;
block.InputPort(1).SamplingMode = 'sample';
block.InputPort(1).SampleTime = [-1 0];

block.OutputPort(1).Dimensions = -1;
block.OutputPort(1).SamplingMode = 'sample';
block.OutputPort(1).SampleTime = [0 0];

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
block_path=gcb;
dot_pos=strfind(block_path,'/');
set_param(block_path(1:dot_pos(length(dot_pos))-1),'MaskDisplay',str);

function CheckParameters(block)
global ext_iqc_setvalues

% Simulink Parameters Check
sta_dea=block.DialogPrm(1).Data;
if ~isa(sta_dea,'double')
    disp_str(9,'start of deadzone','numeric')
end
end_dea=block.DialogPrm(2).Data;
if ~isa(end_dea,'double')
    disp_str(9,'end of deadzone','numeric')
end
if sta_dea>end_dea
    disp_str(39,'"end of deadzone"','<','"start of deadzone"')
end
slo=block.DialogPrm(3).Data;
if ~isa(slo,'double')
    DAStudio.error('slope of deadzone','numeric');
end
if slo<0
    DAStudio.error('slope of deadzone','positive number');
end

% IQC Procedure Parameters Check
block_path=gcb;
ny=strfind(block_path,'/');
block_newpath=[block_path(1:ny-1),'_lin',block_path(ny:ny(end)-1)];
block_name=get_param(gcb,'name');

glparam_name='ext_sim_intenc_deadzone';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','path','newpath'};
alparam={block_name,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);

ext_iqc_setvalues.ext_sim_intenc_deadzone.iqctype{pos}='';

iqc1=block.DialogPrm(4).Data;
iqc2=block.DialogPrm(8).Data;
iqc3=block.DialogPrm(12).Data;
if ~iqc1 && ~iqc2 && ~iqc3
    DAStudio.error('Simulink:block:invalidParameter');
    iqctool('simclear')
end
if iqc1
    iqc1_pole=block.DialogPrm(5).Data;
    if ~isa(iqc1_pole, 'double')
        disp_str(9,'"iqc_dzn_e pole"','numeric')
    end
    iqc1_def_num_h=block.DialogPrm(6).Data;
    if ~isa(iqc1_def_num_h,'double')
        disp_str(9,'"iqc_dzn_e Defines Number ..."','numeric')
    end
    iqc1_out=block.DialogPrm(7).Data;
    for i=1:length(iqc1_out)
        if ~isvarname(iqc1_out{i})
            DAStudio.error('Simulink:block:invalidParameter');
        end
    end
    setvalue(pos,'iqc_dzn_e',iqc1_out,{'pole',iqc1_pole;...
        'def_num_h',iqc1_def_num_h;...
        'slo_bou',slo});
end
if iqc2
    iqc2_pole=block.DialogPrm(9).Data;
    if ~isa(iqc2_pole,'double')
        disp_str(9,'"iqc_dzn_e_odd pole"','numeric')
    end
    iqc2_def_num_h=block.DialogPrm(10).Data;
    if ~isa(iqc2_def_num_h,'double')
        disp_str(9,'"iqc_dzn_e_odd Defines Number ..."','numeric')
    end
    iqc2_out=block.DialogPrm(11).Data;
    for i=1:length(iqc2_out)
        if ~isvarname(iqc2_out{i})
            DAStudio.error('Simulink:block:invalidParameter');
        end
    end
    setvalue(pos,'iqc_dzn_e_odd',iqc2_out,{'pole',iqc2_pole;...
        'def_num_h',iqc2_def_num_h;...
        'slo_bou',slo});
end
if iqc3
    lis_rel=block.DialogPrm(13).Data;
    if ~ischar(lis_rel)
        disp_str(9,'custom IQC relationships','char')
    end
    setvalue(pos,'iqc_custom',{'0'},{'lis_rel',lis_rel});
end

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_intenc_deadzone';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_sim_intenc_deadzone.sampletime{pos}(1:2)=di(1);
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
glparam_name='ext_sim_intenc_deadzone';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_sim_intenc_deadzone.size{pos}(1:2)=di;
% Set compiled dimensions
block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = di;


function Outputs(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_intenc_deadzone';
pos=checkblockpos(block_path,glparam_name);

if ~isfield(ext_iqc_setvalues.ext_sim_intenc_deadzone,'size')
    ext_iqc_setvalues.ext_sim_intenc_deadzone.size{pos}(1:2)=...
        block.InputPort(1).Dimensions;
end

sta_dea=block.DialogPrm(1).Data;
end_dea=block.DialogPrm(2).Data;
slo=block.DialogPrm(3).Data;
ny=size(block.InputPort(1).Data,1);

for i=1:ny,
    if block.InputPort(1).Data(i,1)<end_dea && sta_dea<block.InputPort(1).Data(i,1)
        block.OutputPort(1).Data(i,1)=0;
    elseif end_dea<=block.InputPort(1).Data(i,1)
        block.OutputPort(1).Data(i,1)=slo*(block.InputPort(1).Data(i,1)-end_dea);
    elseif block.InputPort(1).Data(i,1) <= sta_dea
        block.OutputPort(1).Data(i,1)=slo*(block.InputPort(1).Data(i,1)-sta_dea);
    end
end

function Terminate(block)
iqctool('simclear')

function setvalue(pos,iqc_name,out_var,other)
global ext_iqc_setvalues
str=[];
for i=1:length(out_var)
    str=[str,out_var{i},','];
end
try
    ny1=size(ext_iqc_setvalues.ext_sim_intenc_deadzone.iqctype{pos},2);
catch err
    ny1=0;
end
if ~strcmp(str(1:end-1),'0')
    ext_iqc_setvalues.ext_sim_intenc_deadzone.iqctype{pos}{ny1+1}(1,:)={[iqc_name,...
        '.outvariable'],...
        str};
else
    ext_iqc_setvalues.ext_sim_intenc_deadzone.iqctype{pos}{ny1+1}(1,:)={[iqc_name,...
        '.list_iqc_rel'],...
        other{1,2}};
    other=[];
end
for i=1:size(other,1)
    
    ny2=size(ext_iqc_setvalues.ext_sim_intenc_deadzone.iqctype{pos}{ny1+1},1);
    ext_iqc_setvalues.ext_sim_intenc_deadzone.iqctype{pos}{ny1+1}(ny2+1,:)=...
        {[iqc_name,'.',other{i,1}],other{i,2}};
end
