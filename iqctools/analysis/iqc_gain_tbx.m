function gain=iqc_gain_tbx(f,z)
% function gain=iqc_gain_tbx(f,z)
%
% gives best estimate of the L2 gain f -> z in the system of IQC's
% defined using the "abst" environment "iqc"
%
% stores the optimal multipliers in global variable ABSTSOLUTION
% (so that they can be retrieved using "value")
%
% 'f' and 'y' can be of type "inp" or "sgn"
%
% Written by cmj on 2013/4/23


if nargin < 2
    disp_str(4,'Two')
end

global ABST

if ~isfield(ABST,'log'),
    disp_str(12)
end

if ~strcmp(ABST.name,'iqc')
    disp_str(15,'iqc')
end

% symbolic names for interior types:
cst=1;
var=2;
lin=3;
lmi=4;
inp=5;
sgn=6;
vsg=7;
csg=8;
vcs=9;
qfm=10;
iqc=11;
lnk=12;


if (~isa(f,'abst'))||(~isa(z,'abst')),
    disp_str(60)
end
if (ABST.log(double(f),1)~=sgn)&&(ABST.log(double(f),1)~=inp),
    disp_str(61,'First')
end
if (ABST.log(double(z),1)~=sgn)&&(ABST.log(double(z),1)~=inp),
    disp_str(61,'Second')
end

lmitool = ABST.lmitool;

switch lmitool
    case 'lmilab'
        gain = iqc_analysis_l2gain_lmilab(f,z);
    case 'yalmip'
        gain = iqc_analysis_l2gain_yalmip(f,z);
    otherwise
        disp_str(55)
end
