function ext_sim_rate_limiter(block)
setup(block);

function setup(block)

block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 9;

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
str=['plot(0,0,100,100,[18,45,45,18,18],[60,60,26,26,60],',...
    '[62,82,82,62,62],[60,60,26,26,60],[62,0],[44,44],[79,65],',...
    '[44,44],[71,71],[46,57],[68,69,72,74,73,71,69,68,69,72,74],',...
    '[31,29,29,31,34,35,36,38,40,40,39],[31,31],[26,60],[45,38,24],',...
    '[55,55,33],[24,18],[33,33],[99,82],[44,44],[10,12,8,10,10,90,90],',...
    '[44,51,51,44,80,80,44],[7,2],[55,55])'];

block_path=gcb;
dot_pos=strfind(block_path,'/');
set_param(block_path(1:dot_pos(length(dot_pos))-1),'MaskDisplay',str);


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

% IQC Procedure Parameters Check
block_path=gcb;
ny=strfind(block_path,'/');
block_newpath=[block_path(1:ny-1),'_lin',block_path(ny:ny(end)-1)];
block_name=get_param(gcb,'name');

glparam_name='ext_sim_rate_limiter';
pos=checkblockpos(block_path,glparam_name);
alparam_name={'name','path','newpath'};
alparam={block_name,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);

ext_iqc_setvalues.ext_sim_rate_limiter.iqctype{pos}='';

iqc1=block.DialogPrm(5).Data;
iqc2=block.DialogPrm(8).Data;

if ~iqc1 && ~iqc2
    DAStudio.error('Simulink:block:invalidParameter');
    iqctool('simclear')
end
if iqc1
    iqc1_pole=block.DialogPrm(6).Data;
    if ~isa(iqc1_pole,'double')
        disp_str(9,'"iqc_ratelimiter pole"','numeric')
    end
    firstpole=block.DialogPrm(4).Data;
    
    iqc1_out=block.DialogPrm(7).Data;
    for i=1:length(iqc1_out)
        if ~isvarname(iqc1_out{i})
            DAStudio.error('Simulink:block:invalidParameter');
        end
    end
    d_sat=block.DialogPrm(3).Data;
    setvalue(pos,'iqc_ratelimiter',iqc1_out,{'pole',[firstpole iqc1_pole];...
        'k',d_sat});
end

if iqc2
    lis_rel=block.DialogPrm(9).Data;
    if ~ischar(lis_rel)
        disp_str(9,'custom IQC relationships','char')
    end
    setvalue(pos,'iqc_custom',{'0'},{'lis_rel',lis_rel});
end

function SetInputSampleTime(block, idx, di)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_rate_limiter';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_sim_rate_limiter.sampletime{pos}(1:2)=di(1);
if di(1)==0,
    block.InputPort(idx).Sampletime = [0 di(2)];
    block.OutputPort(idx).Sampletime = [0 0];
else
    block.InputPort(idx).Sampletime = di;
    block.OutputPort(idx).Sampletime = [di(1) 0];
end

function SetOutputSampleTime(block, idx, di)


function SetInputPortDims(block, idx, di)
if di~=1
    disp_str(9,'Input of IQC Sim Rate Limiters','scalar');
end
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_rate_limiter';
pos=checkblockpos(block_path,glparam_name);

ext_iqc_setvalues.ext_sim_rate_limiter.size{pos}(1:2)=di;
% Set compiled dimensions
block.InputPort(idx).Dimensions = di;
block.OutputPort(idx).Dimensions = di;


function Outputs(block)
global ext_iqc_setvalues
block_path=gcb;
glparam_name='ext_sim_rate_limiter';
pos=checkblockpos(block_path,glparam_name);

if ~isfield(ext_iqc_setvalues.ext_sim_rate_limiter,'size')
    ext_iqc_setvalues.ext_sim_rate_limiter.size{pos}(1:2)=...
        block.InputPort(1).Dimensions;
end

upp_lim=block.DialogPrm(1).Data;
low_lim=block.DialogPrm(2).Data;
ny=size(block.InputPort(1).Data,1);
d_sat=block.DialogPrm(3).Data;

for i=1:ny,
    if d_sat*block.InputPort(1).Data(i,1)<upp_lim && low_lim<d_sat*block.InputPort(1).Data(i,1)
        block.OutputPort(1).Data(i,1)=d_sat*block.InputPort(1).Data(i,1);
    elseif upp_lim<=d_sat*block.InputPort(1).Data(i,1)
        block.OutputPort(1).Data(i,1)=upp_lim;
    elseif d_sat*block.InputPort(1).Data(i,1) <= low_lim
        block.OutputPort(1).Data(i,1)=low_lim;
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
    ny1=size(ext_iqc_setvalues.ext_sim_rate_limiter.iqctype{pos},2);
catch err
    ny1=0;
end
if ~strcmp(str(1:end-1),'0')
    ext_iqc_setvalues.ext_sim_rate_limiter.iqctype{pos}{ny1+1}(1,:)={[iqc_name,...
        '.outvariable'],...
        str};
else
    ext_iqc_setvalues.ext_sim_rate_limiter.iqctype{pos}{ny1+1}(1,:)={[iqc_name,...
        '.list_iqc_rel'],...
        other{1,2}};
    other=[];
end
for i=1:size(other,1)
    
    ny2=size(ext_iqc_setvalues.ext_sim_rate_limiter.iqctype{pos}{ny1+1},1);
    ext_iqc_setvalues.ext_sim_rate_limiter.iqctype{pos}{ny1+1}(ny2+1,:)=...
        {[iqc_name,'.',other{i,1}],other{i,2}};
end
