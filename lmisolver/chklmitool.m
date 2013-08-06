function chk=chklmitool

% Last modified by cmj on 2013/4/19

% lmitool location
lmilablocation=which('lmivar.m');
yalmiplocation=which('sdpvar.m');

chk=0;

% check
if isempty(yalmiplocation) && isempty(lmilablocation),
    disp_str(44)
elseif ~isempty(yalmiplocation) && ~isempty(lmilablocation),
    chk=3;
elseif ~isempty(lmilablocation),
    chk=1;
elseif ~isempty(yalmiplocation),
    chk=2;
end
