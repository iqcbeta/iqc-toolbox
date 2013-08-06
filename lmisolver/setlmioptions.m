function setlmioptions(varargin)
% function setlmioptions(lmitool,parameter)
% lmitool: 'lmilab' 或是 'yalmip'
% parameter: 1. 當 lmitool 為 'lmilab' 時，parameter 是 1x5 數值矩陣
%            2. 當 lmitool 為 'yalmip' 時，parameter 如同 sdpsettings 的設定
%
% Ex. 1. 使用 lmilab 且使用預設值
%        setlmioptions('lmilab')
% Ex. 2. 使用 yalmip 且由工具自動選擇解法器
%        setlimoptions('yalmip')
% Ex. 3. 使用 lmilab 且使用設定 [0 300 0 0 777]
%        setlmioptions('lmilab',[0 300 0 0 777])
% Ex. 4. 使用 yalmip 且採用 sedumi 解法器，並且在設定 sedumi.maxit 值為 300
%        setlmioptions('yalmip','solver','sedumi','sedumi.maxit',300)
%
% 設定錯誤將自動忽略本次設定並產生警告
% 祥細參考 yalmip 的 sdpsettings 設定
%
% Last modified by cmj on 2013/4/22

% 輸入參數判斷
if nargin < 1
    disp_str(4,'One')
elseif nargin == 1
    lmioptions='';
end
lmitool=varargin{1};

global ABST

% 環境判斷
if isempty(ABST)
    disp_str(11)
    return
end

% lmitool 檢查及輸入參數檢查
chklmitoolvalue=chklmitool;
if chklmitoolvalue==0 % 沒有可用的 LMI 工具
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
            % 參數輸入錯誤
            disp_str(46)
            return
        end
        lmioptions = varargin{2};
    end
elseif strcmpi(lmitool,'yalmip') &&...
        (chklmitoolvalue==3 || chklmitoolvalue==2)
    uselmitool = 'yalmip';
    
    % parameter 個數非偶數，且第一參數名不是字元
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

% 參數設定顯示
switch uselmitool
    case 'lmilab'
        if isempty(lmioptions) % 不設定 lmi 參數，使用預設 [0 200 0 0 0]
            lmioptions = [0 200 0 0 0];
            disp_str(48,'lmilab')
            displaylmilab(lmioptions)
        else
            displaylmilab(lmioptions) % 使用者設定
        end
    case 'yalmip'
        if isempty(lmioptions) % 不設定 yalmip 參數，使用 yalmip 預設
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
