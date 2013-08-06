function IQCSimSystem(block)

setup(block);

function setup(block)
block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 5;

sam_tim=block.DialogPrm(5).Data;
A=block.DialogPrm(1).Data;
C=block.DialogPrm(3).Data;

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions = -1;
block.InputPort(1).DirectFeedthrough = true;

block.OutputPort(1).Dimensions = size(C,1);

block.SampleTimes = [sam_tim 0];
if sam_tim~=0
    block.RegBlockMethod('Update',@Update);
    block.RegBlockMethod('PostPropagationSetup',@DoPostPropSetup);
else
    block.NumContStates = size(A,1);
    block.RegBlockMethod('Derivatives',@Derivative);
end

block.RegBlockMethod('InitializeConditions',@InitConditions);
block.RegBlockMethod('Outputs',@Output);

function DoPostPropSetup(block)
A=block.DialogPrm(1).Data;
%% Setup Dwork
block.NumDworks = 1;
block.Dwork(1).Name = 'X';
block.Dwork(1).Dimensions      = size(A,1);
block.Dwork(1).DatatypeID      = 0;
block.Dwork(1).Complexity      = 'Real';
block.Dwork(1).UsedAsDiscState = true;

function InitConditions(block)
sam_tim=block.DialogPrm(5).Data;
A=block.DialogPrm(1).Data;

if sam_tim~=0
    block.Dwork(1).Data = zeros(size(A,1),1);
else
    block.ContStates.Data = zeros(size(A,1),1);
end

function Derivative(block)
A=block.DialogPrm(1).Data;
B=block.DialogPrm(2).Data;
u =  block.InputPort(1).Data;

block.Derivatives.Data = A*block.ContStates.Data+B*u;

function Update(block)
A=block.DialogPrm(1).Data;
B=block.DialogPrm(2).Data;
u =  block.InputPort(1).Data;

block.Dwork(1).Data = A*block.Dwork(1).Data+B*u;

function Output(block)
C=block.DialogPrm(3).Data;
D=block.DialogPrm(4).Data;
sam_tim=block.DialogPrm(5).Data;
u =  block.InputPort(1).Data;

if sam_tim~=0
    block.OutputPort(1).Data = C*block.Dwork(1).Data+D*u;
else
    block.OutputPort(1).Data = C*block.ContStates.Data+D*u;
end
