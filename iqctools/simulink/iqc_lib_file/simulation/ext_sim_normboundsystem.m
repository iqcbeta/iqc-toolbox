function ext_sim_normboundsystem(block)
setup(block);

function setup(block)

block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 18;

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions = -1;
block.InputPort(1).DirectFeedthrough = true;
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
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('InitializeConditions',@InitializeConditions);

function InitializeConditions(block)
str=['disp(''\fontsize{14}||\Delta|| < k'', ''texmode'',''on'');'];

block_path=gcb;
dot_pos=strfind(block_path,'/');
set_param(block_path(1:dot_pos(length(dot_pos))-1),'MaskDisplay',str);


function CheckParameters(block)
global ext_iqc_setvalues

% Simulink Parameters Check
upp_bou=block.DialogPrm(5).Data;

% IQC Procedure Parameters Check
block_path=gcb;
ny=strfind(block_path,'/');
block_newpath=[block_path(1:ny-1),'_lin',block_path(ny:ny(end)-1)];
block_name=get_param(gcb,'name');

glparam_name='ext_sim_normboundsystem';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','path','newpath'};
alparam={block_name,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);

ext_iqc_setvalues.ext_sim_normboundsystem.iqctype{pos}='';

iqc1=block.DialogPrm(9).Data;
iqc2=block.DialogPrm(12).Data;
if ~iqc1 && ~iqc2
    DAStudio.error('Simulink:block:invalidParameter');
    iqctool('simclear')
end
if iqc1
    iqc1_pole=block.DialogPrm(10).Data;
    if ~isa(iqc1_pole,'double')
        disp_str(9,'"iqc_ltiunmod pole"','numeric')
    end
    iqc1_out=block.DialogPrm(11).Data;
    for i=1:length(iqc1_out)
        if ~isvarname(iqc1_out{i})
            DAStudio.error('Simulink:block:invalidParameter');
        end
    end
    num_out=block.DialogPrm(3).Data;
    setvalue(pos,'iqc_ltiunmod',iqc1_out,{'pole',iqc1_pole;...
        'upp_bou',upp_bou;...
        'out_dim',num_out});
end
if iqc2
    lis_rel=block.DialogPrm(13).Data;
    if ~ischar(lis_rel)
        disp_str(9,'custom IQC relationships','char')
    end
    setvalue(pos,'iqc_custom',{'0'},{'lis_rel',lis_rel});
end

function Start(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_normboundsystem';
pos=checkblockpos(block_path,glparam_name);

sam_tim=block.DialogPrm(6).Data;
A=block.DialogPrm(15).Data;
B=block.DialogPrm(16).Data;
C=block.DialogPrm(17).Data;
D=block.DialogPrm(18).Data;
if sam_tim==0
    G=ss(ss(A,B,C,D));
else
    G=ss(ss(A,B,C,D,sam_tim));
end
block_name=gcb;
dot_pos=strfind(gcb,'/');
n=length(dot_pos);
block_name=block_name(dot_pos(n-1)+1:dot_pos(n)-1);
fprintf('The "%s" Simulink Parameter\n',block_name)
fprintf('The Random System:\n');
tf(G)
ext_iqc_setvalues.ext_sim_normboundsystem.system{pos}=G;
var_nam=block.DialogPrm(8).Data;
if ~isempty(var_nam)
    eval([var_nam,'=G;']);
    eval(['global ',var_nam,';']);
end


function SetInputSampleTime(block, idx, di)
sam_tim=block.DialogPrm(6).Data;
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_normboundsystem';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_sim_normboundsystem.sampletime{pos}(1:2)=di(1);

if sam_tim==0,
    block.InputPort(idx).Sampletime = [0 di(2)];
    block.OutputPort(idx).Sampletime = [0 0];
else
    block.InputPort(idx).Sampletime = di;
    %     block.OutputPort(idx).Sampletime = [sam_tim/100 di(2)];
    block.OutputPort(idx).Sampletime = [sam_tim di(2)];
end

function SetOutputSampleTime(block, idx, di)


function SetInputPortDims(block, idx, di)
num_out=block.DialogPrm(3).Data;
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_normboundsystem';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_sim_normboundsystem.size{pos}(1)=di;
ext_iqc_setvalues.ext_sim_normboundsystem.size{pos}(2)=num_out;
% Set compiled dimensions
block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = num_out;


function Outputs(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_normboundsystem';
pos=checkblockpos(block_path,glparam_name);

num_out=block.DialogPrm(3).Data;

if ~isfield(ext_iqc_setvalues.ext_sim_normboundsystem,'size')
    ext_iqc_setvalues.ext_sim_normboundsystem.size{pos}(1)=...
        block.InputPort(1).Dimensions;
    ext_iqc_setvalues.ext_sim_normboundsystem.size{pos}(2)=num_out;
end
block.OutputPort(1).Data=zeros(block.OutputPort(1).Dimensions,1);
%endfunction

function Terminate(block)
iqctool('simclear')
%end Terminate

function setvalue(pos,iqc_name,out_var,other)
global ext_iqc_setvalues
str=[];
for i=1:length(out_var)
    str=[str,out_var{i},','];
end
try
    ny1=size(ext_iqc_setvalues.ext_sim_normboundsystem.iqctype{pos},2);
catch err
    ny1=0;
end
if ~strcmp(str(1:end-1),'0')
    ext_iqc_setvalues.ext_sim_normboundsystem.iqctype{pos}{ny1+1}(1,:)={[iqc_name,...
        '.outvariable'],...
        str};
else
    ext_iqc_setvalues.ext_sim_normboundsystem.iqctype{pos}{ny1+1}(1,:)={[iqc_name,...
        '.list_iqc_rel'],...
        other{1,2}};
    other=[];
end
for i=1:size(other,1)
    
    ny2=size(ext_iqc_setvalues.ext_sim_normboundsystem.iqctype{pos}{ny1+1},1);
    ext_iqc_setvalues.ext_sim_normboundsystem.iqctype{pos}{ny1+1}(ny2+1,:)=...
        {[iqc_name,'.',other{i,1}],other{i,2}};
end
