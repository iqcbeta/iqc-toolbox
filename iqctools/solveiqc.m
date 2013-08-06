function varargout = solveiqc(problem,method,parameter)

if nargin == 0
    disp_str(1)
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

switch problem
    case 'analysis'
        switch method
            case 'l2gain'
                if length(parameter)~=2 || ~iscell(parameter)
                    disp_str(69,'parameter','{signal1,signal2}')
                end
                f=parameter{1};
                z=parameter{2};
                
                %% ¦P iqc_gain_tbx
                if (~isa(f,'abst'))||(~isa(z,'abst')),
                    disp_str(60)
                end
                if (ABST.log(double(f),1)~=sgn)&&...
                        (ABST.log(double(f),1)~=inp),
                    disp_str(61,'First')
                end
                if (ABST.log(double(z),1)~=sgn)&&...
                        (ABST.log(double(z),1)~=inp),
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
                %%
                varargout = {gain};
                return
            otherwise
                disp_str(69,'method','''l2gain''')
        end
    case 'estimation'
        switch method
            case 'l2gain-rep'
                [gain,estimation]=setest('l2gain-rep',parameter);
            case 'l2gain-eli'
                [gain,estimation]=setest('l2gain-eli',parameter);
            otherwise
                disp_str(69,'method','''l2gain-rep'',''l2gain-eli''')
        end
        varargout = {gain,estimation};
        return
    otherwise
        disp_str(69,'problem','''analysis'', ''estimation''')
end

function [gain,estimation]=setest(method,parameter)

global ABST

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

if length(parameter) < 6
    parameter{6}=1;
end
if length(parameter) < 5
    disp_str(69,'parameter','{w,v,y,p,q,W}')
end
w=parameter{1};
v=parameter{2};
y=parameter{3};
p=parameter{4};
q=parameter{5};
W=parameter{6};
if (~isa(w,'abst'))||(~isa(v,'abst'))||(~isa(y,'abst'))||...
        (~isa(p,'abst'))||(~isa(q,'abst'))
    disp_str(60)
end

if (ABST.log(double(w),1)~=sgn)&&...
        (ABST.log(double(w),1)~=inp),
    disp_str(61,'First')
end
if (ABST.log(double(v),1)~=sgn)&&...
        (ABST.log(double(v),1)~=inp),
    disp_str(61,'Second')
end
if (ABST.log(double(y),1)~=sgn)&&...
        (ABST.log(double(y),1)~=inp),
    disp_str(61,'Third')
end
if (ABST.log(double(p),1)~=sgn)&&...
        (ABST.log(double(p),1)~=inp),
    disp_str(61,'Fourth')
end
if (ABST.log(double(q),1)~=sgn)&&...
        (ABST.log(double(q),1)~=inp),
    disp_str(61,'Fifth')
end

lmitool = ABST.lmitool;
systemtype = ABST.systemtype;

chk_par = [systemtype,'_',lmitool,'_',method];

switch chk_par
    case 'continuous_lmilab_l2gain-rep'
        [gain,estimation]=iqc_estimation_l2gain_CTrep_lmilab(w,v,y,p,q,W);
    case 'continuous_yalmip_l2gain-rep'
        [gain,estimation]=iqc_estimation_l2gain_CTrep_yalmip(w,v,y,p,q,W);
    case 'continuous_lmilab_l2gain-eli'
        [gain,estimation]=iqc_estimation_l2gain_CTeli_lmilab(w,v,y,p,q,W);
    case 'continuous_yalmip_l2gain-eli'
        [gain,estimation]=iqc_estimation_l2gain_CTeli_yalmip(w,v,y,p,q,W);
        
    case 'discrete_lmilab_l2gain-rep'
        [gain,estimation]=iqc_estimation_l2gain_DTrep_lmilab(w,v,y,p,q,W);
    case 'discrete_yalmip_l2gain-rep'
        [gain,estimation]=iqc_estimation_l2gain_DTrep_yalmip(w,v,y,p,q,W);
    case 'discrete_lmilab_l2gain-eli'
        [gain,estimation]=iqc_estimation_l2gain_DTeli_lmilab(w,v,y,p,q,W);
    case 'discrete_yalmip_l2gain-eli'
        [gain,estimation]=iqc_estimation_l2gain_DTeli_yalmip(w,v,y,p,q,W);
end
