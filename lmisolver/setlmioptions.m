function setlmioptions(varargin)
% function setlmioptions(lmitool,parameter)
% lmitool: 'lmilab' �άO 'yalmip'
% parameter: 1. �� lmitool �� 'lmilab' �ɡAparameter �O 1x5 �ƭȯx�}
%            2. �� lmitool �� 'yalmip' �ɡAparameter �p�P sdpsettings ���]�w
%
% Ex. 1. �ϥ� lmilab �B�ϥιw�]��
%        setlmioptions('lmilab')
% Ex. 2. �ϥ� yalmip �B�Ѥu��۰ʿ�ܸѪk��
%        setlimoptions('yalmip')
% Ex. 3. �ϥ� lmilab �B�ϥγ]�w [0 300 0 0 777]
%        setlmioptions('lmilab',[0 300 0 0 777])
% Ex. 4. �ϥ� yalmip �B�ĥ� sedumi �Ѫk���A�åB�b�]�w sedumi.maxit �Ȭ� 300
%        setlmioptions('yalmip','solver','sedumi','sedumi.maxit',300)
%
% �]�w���~�N�۰ʩ��������]�w�ò���ĵ�i
% ���ӰѦ� yalmip �� sdpsettings �]�w
%
% Last modified by cmj on 2013/4/22

% ��J�ѼƧP�_
if nargin < 1
    disp_str(4,'One')
elseif nargin == 1
    lmioptions='';
end
lmitool=varargin{1};

global ABST

% ���ҧP�_
if isempty(ABST)
    disp_str(11)
    return
end

% lmitool �ˬd�ο�J�Ѽ��ˬd
chklmitoolvalue=chklmitool;
if chklmitoolvalue==0 % �S���i�Ϊ� LMI �u��
    return
end
if strcmpi(lmitool,'lmilab') && (chklmitoolvalue==3 || chklmitoolvalue==1)
    uselmitool = 'lmilab';
    if nargin > 2
        disp_str(2)
        return
    end
    if nargin==2
        if ~isa(varargin{2},'double') ||...
                (size(varargin{2},1)~=1 || size(varargin{2},2)~=5)
            % �Ѽƿ�J���~
            disp_str(46)
            return
        end
        lmioptions = varargin{2};
    end
elseif strcmpi(lmitool,'yalmip') &&...
        (chklmitoolvalue==3 || chklmitoolvalue==2)
    uselmitool = 'yalmip';
    
    % parameter �ӼƫD���ơA�B�Ĥ@�ѼƦW���O�r��
    for i1=2:2:nargin
        if ~isa(varargin{i1},'char')
            disp_str(47)
            return
        end
    end
    if rem(length(varargin)-1,2)==1
        disp_str(47)
        return
    end
    lmioptions = varargin(2:end);
else
    disp_str(45,lmitool)
    return
end

% �ѼƳ]�w���
switch uselmitool
    case 'lmilab'
        if isempty(lmioptions) % ���]�w lmi �ѼơA�ϥιw�] [0 200 0 0 0]
            lmioptions = [0 200 0 0 0];
            disp_str(48,'lmilab')
            displaylmilab(lmioptions)
        else
            displaylmilab(lmioptions) % �ϥΪ̳]�w
        end
    case 'yalmip'
        if isempty(lmioptions) % ���]�w yalmip �ѼơA�ϥ� yalmip �w�]
            disp_str(48,'yalmip')
        else
            options_temp='';
            for i1=1:length(lmioptions)
                if isa(lmioptions{i1},'char')
                    options_temp=[options_temp,'''',lmioptions{i1},''',']; %#ok<*AGROW>
                else
                    options_temp=[options_temp,...
                        mat2str(lmioptions{i1}),','];
                end
            end
            displayyalmip(lmioptions)
            lmioptions=options_temp(1:end-1);
        end
end
ABST.lmitool = uselmitool;
ABST.lmiparameter = lmioptions;


%%
function displaylmilab(lmioptions)
%%
s=[];
%     12345678901234567890123456789012345678901234567890
s='   LMI Control Toolbox      options:              ';
s=[s;'**************************************************'];
s=[s;' 1. Relative accuracy:        ' sft(lmioptions(1),20)];
s=[s;' 2. Maximal # of iterations:  ' sft(lmioptions(2),20)];
s=[s;' 3. Feasibility radius:       ' sft(lmioptions(3),20)];
s=[s;' 4. Slow progress #:          ' sft(lmioptions(4),20)];
s=[s;' 5. Turn-off trace:           ' sft(lmioptions(5),20)];
s=[s;'**************************************************'];
disp(s)

%%
function displayyalmip(lmioptions)
%%
s=[];
ns=1;
%     12345678901234567890123456789012345678901234567890
s='   YALMIP Toolbox           options:              ';
s=[s;'**************************************************'];
for i1=1:2:length(lmioptions)
    s=[s;[sft(ns,2),'. ',sft(lmioptions{i1},26),sft(lmioptions{i1+1},20)]];
    ns=ns+1;
end
s=[s;'**************************************************'];
disp(s)
