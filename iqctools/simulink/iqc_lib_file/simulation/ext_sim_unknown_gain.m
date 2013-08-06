function ext_sim_unknown_gain(block)
setup(block);

function setup(block)

block.NumInputPorts  = 1;
block.NumOutputPorts = 2;
block.NumDialogPrms  = 13;

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions = -1;
block.InputPort(1).DirectFeedthrough = true;
block.InputPort(1).SamplingMode = 'sample';
block.InputPort(1).SampleTime = [-1 0];

block.OutputPort(1).Dimensions = -1;
block.OutputPort(1).SamplingMode = 'sample';
block.OutputPort(1).SampleTime = [-1 0];

block.OutputPort(2).Dimensions = -1;
block.OutputPort(2).SamplingMode = 'sample';
block.OutputPort(2).SampleTime = [0 0];

block.RegBlockMethod('SetInputPortDimensions',@SetInputPortDims);
block.RegBlockMethod('SetInputPortSampleTime',@SetInputSampleTime);
block.RegBlockMethod('SetOutputPortSampleTime',@SetOutputSampleTime);
block.RegBlockMethod('CheckParameters',@CheckParameters);
block.RegBlockMethod('Outputs', @Outputs);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Terminate', @Terminate);
block.RegBlockMethod('InitializeConditions',@InitializeConditions);

function InitializeConditions(block)
str=['disp(''\fontsize{14}|d| < 1'', ''texmode'',''on'');',...
    'port_label(''input'',1,''in'');',...
    'port_label(''output'',1,''out'');',...
    'port_label(''output'',2,''|d|'');'];
set_param(gcb,'MaskDisplay',str);

function CheckParameters(block)
global ext_iqc_setvalues

% IQC Procedure Parameters Check
block_path=gcb;
ny=strfind(block_path,'/');
block_newpath=[block_path(1:ny-1),'_lin',block_path(ny:size(block_path,2))];
block_name=get_param(gcb,'name');

glparam_name='ext_sim_unknown_gain';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','path','newpath'};
alparam={block_name,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);

ext_iqc_setvalues.ext_sim_unknown_gain.iqctype{pos}='';

iqc1=block.DialogPrm(4).Data;
iqc2=block.DialogPrm(7).Data;
iqc3=block.DialogPrm(9).Data;
iqc4=block.DialogPrm(11).Data;
if ~iqc1 && ~iqc2 && ~iqc3 && ~iqc4
    DAStudio.error('Simulink:block:invalidParameter');
    iqctool('simclear')
end
if iqc1
    iqc1_pole=block.DialogPrm(5).Data;
    if ~isa(iqc1_pole,'double')
        disp_str(9,'"iqc_ltigain pole"','numeric')
    end
    iqc1_out=block.DialogPrm(6).Data;
    for i=1:length(iqc1_out)
        if ~isvarname(iqc1_out{i})
            DAStudio.error('Simulink:block:invalidParameter');
        end
    end
    setvalue(pos,'iqc_ltigain',iqc1_out,{'pole',iqc1_pole});
end
if iqc2
    iqc2_out=block.DialogPrm(8).Data;
    for i=1:length(iqc2_out)
        if ~isvarname(iqc2_out{i})
            DAStudio.error('Simulink:block:invalidParameter');
        end
    end
    setvalue(pos,'iqc_ltvnorm',iqc2_out,{'out_dim',-1;...
        'max_amp',1});
end
if iqc3
    iqc3_out=block.DialogPrm(10).Data;
    for i=1:length(iqc3_out)
        if ~isvarname(iqc3_out{i})
            DAStudio.error('Simulink:block:invalidParameter');
        end
    end
    setvalue(pos,'iqc_tvscalar',iqc3_out,{'upp_bou',1});
end
if iqc4
    lis_rel=block.DialogPrm(12).Data;
    if ~ischar(lis_rel)
        disp_str(9,'custom IQC relationships','char')
    end
    setvalue(pos,'iqc_custom',{'0'},{'lis_rel',lis_rel});
end

function Start(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_unknown_gain';
pos=checkblockpos(block_path,glparam_name);

rand_cons=block.DialogPrm(13).Data;

ext_iqc_setvalues.ext_sim_unknown_gain.rand_cons{pos}=rand_cons;
inputdim=block.InputPort(1).Dimensions;
block_name=gcb;
dot_pos=strfind(gcb,'/');
n=length(dot_pos);
block_name=block_name(dot_pos(n)+1:length(block_name));
fprintf('The "%s" Simulink Parameter\n',block_name)
if inputdim==size(rand_cons,2)
    for i=1:size(rand_cons,2)
        fprintf('The #%d Random Constant: %2.3f \n',i,rand_cons(1,i));
    end
elseif inputdim~=1 && size(rand_cons,2)==1
    for i=1:inputdim
        fprintf('The #%d Random Constant: %2.3f \n',...
            i,rand_cons(1,1));
    end
else
    error('The Number of Random ~= Input dimension')
end

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_unknown_gain';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_sim_unknown_gain.sampletime{pos}(1:2)=di(1);
if di(1)==0,
    block.InputPort(idx).Sampletime = [0 di(2)];
    block.OutputPort(idx).Sampletime = [0 0];
else
    block.InputPort(idx).Sampletime = di;
    block.OutputPort(idx).Sampletime = [di(1) 0];
end
ext_iqc_setvalues.ext_sim_unknown_gain.size{pos}(1:2)=...
    block.InputPort(1).Dimensions;
function SetOutputSampleTime(block, idx, di)


function SetInputPortDims(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_unknown_gain';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_sim_unknown_gain.size{pos}(1:2)=di;
iqc2=block.DialogPrm(7).Data;
if iqc2
    maxi1=length(ext_iqc_setvalues.ext_sim_unknown_gain.iqctype{pos});
    for i1=1:maxi1
        maxi2=size(ext_iqc_setvalues.ext_sim_unknown_gain.iqctype{pos}{1,i1},1);
        for i2=1:maxi2
            if strcmp(ext_iqc_setvalues.ext_sim_unknown_gain.iqctype{pos}{1,i1}{i2,1},...
                    'iqc_ltvnorm.out_dim')
                ext_iqc_setvalues.ext_sim_unknown_gain.iqctype{pos}{1,i1}{i2,2}=di;
                break;
            end
        end
    end
end
rand_cons=block.DialogPrm(13).Data;


block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = di;
    block.OutputPort(2).Dimensions = length(rand_cons);

    
function Outputs(block)
block_path=gcb;
glparam_name='ext_sim_unknown_gain';
pos=checkblockpos(block_path,glparam_name);

global ext_iqc_setvalues
if ~isfield(ext_iqc_setvalues.ext_sim_unknown_gain,'size')
    ext_iqc_setvalues.ext_sim_unknown_gain.size{pos}(1:2)=...
        block.InputPort(1).Dimensions;
end
rand_cons=ext_iqc_setvalues.ext_sim_unknown_gain.rand_cons{pos};
ny=length(rand_cons);
% rand_cons=reshape(rand_cons,ny,1);
% outmat=zeros(block.OutputPort(1).Dimensions,block.InputPort(1).Dimensions);
% for i1=1:block.InputPort(1).Dimensions
%     outmat((i1-1)*ny+1:i1*ny,i1)=rand_cons;
% end
% block.OutputPort(1).Data=outmat*block.InputPort(1).Data;

if ny==1 && block.OutputPort(1).Dimensions~=1
    for i=1:block.OutputPort(1).Dimensions
        block.OutputPort(1).Data(i,1)=rand_cons*block.InputPort(1).Data(i,1);
    end
elseif ny==block.OutputPort(1).Dimensions
    for i=1:block.OutputPort(1).Dimensions
        block.OutputPort(1).Data(i,1)=rand_cons(1,i)*block.InputPort(1).Data(i,1);
    end
end
block.OutputPort(2).Data=rand_cons';

function Terminate(block)
iqctool('simclear')


function setvalue(pos,iqc_name,out_var,other)
global ext_iqc_setvalues
str=[];
for i=1:length(out_var)
    str=[str,out_var{i},','];
end
try
    ny1=size(ext_iqc_setvalues.ext_sim_unknown_gain.iqctype{pos},2);
catch err
    ny1=0;
end
if ~strcmp(str(1:end-1),'0')
    ext_iqc_setvalues.ext_sim_unknown_gain.iqctype{pos}{ny1+1}(1,:)={[iqc_name,...
        '.outvariable'],...
        str};
else
    ext_iqc_setvalues.ext_sim_unknown_gain.iqctype{pos}{ny1+1}(1,:)={[iqc_name,...
        '.list_iqc_rel'],...
        other{1,2}};
    other=[];
end
for i=1:size(other,1)
    ny2=size(ext_iqc_setvalues.ext_sim_unknown_gain.iqctype{pos}{ny1+1},1);
    ext_iqc_setvalues.ext_sim_unknown_gain.iqctype{pos}{ny1+1}(ny2+1,:)=...
        {[iqc_name,'.',other{i,1}],other{i,2}};
end
