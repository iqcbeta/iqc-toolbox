function ext_iqc_performance(block)
setup(block);

function setup(block)

block.NumInputPorts  = 1;
block.NumOutputPorts = 0;
block.NumDialogPrms  = 6;

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions = -1;
block.InputPort(1).DirectFeedthrough = false;
block.InputPort(1).SamplingMode = 'sample';
block.InputPort(1).SampleTime = [-1 0];

block.RegBlockMethod('SetInputPortDimensions',@SetInputPortDims);
block.RegBlockMethod('Outputs', @Outputs);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Terminate', @Terminate);
block.RegBlockMethod('SetInputPortSampleTime',@SetInputSampleTime);
block.RegBlockMethod('InitializeConditions',@InitializeConditions);

function InitializeConditions(block)
str='disp(''In -> Out \n performance \n L2 gain:\n Unknown'');';
set_param(gcb,'MaskDisplay',str);

function Start(block)
global ext_iqc_setvalues
save_lin_com=block.DialogPrm(2).Data;
uselmitool=block.DialogPrm(5).Data;
lmi_par=block.DialogPrm(6).Data;

ext_iqc_setvalues.runguiiqc=0;

block_path=gcb;
ny=strfind(block_path,'/');
block_newpath=[block_path(1:ny-1),'_lin',block_path(ny:size(block_path,2))];

block_name=get_param(gcb,'name');

glparam_name='ext_iqc_performance';
pos=checkblockpos(block_path,glparam_name);
ext_iqc_setvalues.ext_iqc_performance.jump_idx=0;

if pos>1
    disp_str(73)
    iqctool('simclear');
end

alparam_name={'name','save_lin_com',...
    'uselmitool','lmi_par','path','newpath'};
alparam={block_name,save_lin_com,...
    uselmitool,lmi_par,block_path,block_newpath};
setblockparam(pos,glparam_name,alparam_name,alparam);
ext_iqc_setvalues.ext_iqc_performance.size{pos}(1)=...
    block.InputPort(1).Dimensions;
ext_iqc_setvalues.ext_iqc_performance.sampletime{pos}(1)=...
    block.InputPort(1).Sampletime(1);

function SetInputSampleTime(block, idx, di)
block.InputPort(idx).Sampletime = di;

function SetInputPortDims(block, idx, di)
block.InputPort(idx).Dimensions = di;

function Outputs(block)
global ABSTSOLUTION
global ext_iqc_setvalues

imp_met=block.DialogPrm(1).Data;
if ext_iqc_setvalues.runguiiqc==0
    ext_iqc_setvalues.runguiiqc=1;
    sys_name=gcs;
    switch imp_met
        case 3
            str='disp(''In -> Out \n performance \n L2 gain:\n Unknown'');';
            set_param(gcb,'MaskDisplay',str);
            stopvalue=findstopsimblk;
        otherwise
            isiqcbode=block.DialogPrm(3).Data;
            fre_ran=block.DialogPrm(4).Data;
            gain=iqc_gui(gcs);
            if isiqcbode
                switch length(fre_ran)
                    case 0
                        isinput=0;
                    case 1
                        isinput=1;
                        input=num2str(fre_ran);
                    case 2
                        isinput=1;
                        input=[num2str(fre_ran(1)),',',num2str(fre_ran(2))];
                end
                if isempty(ABSTSOLUTION)
                    iqc_value
                end
                if ~isinput
                    iqc_bode
                else
                    eval(['iqc_bode(',input,')'])
                end
            end
            block_path=ext_iqc_setvalues.ext_iqc_performance.path{1};
            if ~isempty(gain)
                str=['fprintf(''In -> Out \n performance \n L2 gain:\n%4.4f'',',...
                    num2str(gain),');'];
            else
                str='disp(''In -> Out \n performance \n L2 gain:\n fail'');';
            end
            set_param(block_path,'MaskDisplay',str);
            if imp_met==2
                set_param(sys_name,'SimulationCommand', 'stop')
                stopvalue=0;
            else
                stopvalue=findstopsimblk;
            end
    end
    
    save_system(sys_name,[],'OverwriteIfChangedOnDisk',true);
    if stopvalue % 由於包含純計算方塊(Procedure Block)所以停止模擬
        warning('on')
        warning('Stop simulation!')
        warning('off')
        set_param(sys_name,'SimulationCommand', 'stop')
    end
end

function Terminate(block)
iqctool('simclear')


function out_var=findstopsimblk
m_dir=which('ext_iqc_cdelay.m');
l=dir(m_dir(1:length(m_dir)-17));
l(1:2,:)='';
global ext_iqc_setvalues

ext=fieldnames(ext_iqc_setvalues);

for i1=1:length(ext)
    for i2=1:length(l)
        if strcmp([ext{i1,1},'.m'],eval('l(i2,1).name'))
            out_var=1;
            return
        end
    end
end
out_var=0;