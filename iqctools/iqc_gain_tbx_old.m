function gain=lmi_gain_tbx_old(f,z)
% function gain=lmi_gain_tbx_old(f,z)
%
% gives best estimate of the L2 gain f->z in the system of IQC's
% defined using the "abst" environment "iqc"
%
% stores the optimal multipliers in global variable ABSTSOLUTION
% (so that they can be retrieved using "value")
%
% f,y  can be of type "inp" or "sgn" 
%
% this function uses the LMI Control Toolbox by The MathWorks, Inc.
% if ABST.options(5)>100, the resulting LMI script is
% exported in the format of two files: lmi_exe.m and lmi_exe.mat
% (containing the script and the datafile respectively)
%
% Written by ameg@mit.edu, last modified November 10, 1997

% first, process the iqc abst log
E=iqc_extract_old(f,z);


% vrb=0 to suppress the output
vrb=(E.options(5)==0)|(E.options(5)==777);

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

if nargin~=2, error('two inputs required'); end
if (~isa(f,'abst'))|(~isa(z,'abst')),
   error('all inputs must be from the "abst" class')
end
if (E.T(double(f))~=sgn)&(E.T(double(f))~=inp),
   error('First argument must be input or signal')
end
if (E.T(double(z))~=sgn)&(E.T(double(z))~=inp),
   error('Second argument must be input or signal')
end
if vrb,
   if E.options(5)>100,
      disp('iqc_gain_tbx  (exporting the LMI script to lmi_exe.m) ...')
   else
      disp('iqc_gain_tbx ...')
   end
end

export_date=date;
str{1}=['%% Created by iqc_gain_tbx on ' export_date];
str{2}='\n%% Load coefficient matrices ...';
str{3}=['load lmi_exe'];

% initialize the LMI Control Toolbox
str{4}='\n%% Initialize the LMI Lab';
str{5}='setlmis([]);';   
eval(str{5})

% now define all variables
if vrb,
   disp('  defining the original variables ...')
end
kvar=find(E.T==var);
str{6}='\n%% Define multiplier variables ...';
sc=7;              % string counter
for k=1:length(kvar),
   ks=num2str(k);
   eval(['struct' ks '=E.L{' num2str(kvar(k)) '};'])
   str{sc}=['x' ks '=lmivar(3,struct' ks ');'];
   eval(str{sc});
   sc=sc+1;
end
      
% now define all LMI's
if vrb,
   disp('  defining the non-KYP LMIs ...')
end
kvar=find(E.T==lmi);
str{sc}='\n%% Define non-KYP lmis ...';
sc=sc+1;
for k=1:length(kvar),
   kk=kvar(k);
   [nn,mm]=size(E.AB{kk});    % dimensions
   ks=num2str(k);
   nns=num2str(nn);
   vs=mm-nn;                       % LMI size
   vss=num2str(vs);
   if nn>0,
      eval(['state' ks '=E.AB{kk};']);
      str{sc}=['p' ks '=lmivar(1,[' nns ' 1]);'];
      eval(str{sc});
      sc=sc+1;
      str{sc}=['lmiterm([' ks ' 1 1 p' ks ...
         '],[eye(' nns ');zeros(' vss ',' nns ')],state' ks ...
         ',''s'');'];
      eval(str{sc});
      sc=sc+1;
   end
   nh=0;
   ng=sum(E.X{kk}(3,:));
   for i=1:size(E.X{kk},2),
      is=num2str(i);
      nvs=num2str(E.X{kk}(1,i));   % variable number
      eval(['l_' ks '_' is '=E.C{kk}(ng+1:ng+E.X{kk}(2,i),:)'';'])
      eval(['r_' ks '_' is '=E.C{kk}(nh+1:nh+E.X{kk}(3,i),:);'])
      str{sc}=['lmiterm([' ks ' 1 1 ' nvs '],l_' ks ...
            '_' is ',r_' ks '_' is ',''s'');'];
      eval(str{sc});
      sc=sc+1;
      nh=nh+E.X{kk}(3,i);
      ng=ng+E.X{kk}(2,i);
   end
end

% now we define the "main" lmi 
if vrb,
   disp('  defining the KYP LMIs ...')
end
str{sc}='\n%% KYP LMI terms...';
sc=sc+1;
nm=E.nlmi+1;
nms=num2str(nm);
mm=E.ninputs;
mms=num2str(mm);
if E.nstates>0,
   nss=num2str(E.nstates);
   str{sc}=['P=lmivar(1,[' nss ' 1]);'];
   eval(str{sc});
   sc=sc+1;
   state0=E.ab;
   % lmiterm([-nm-1 1 1 P],1,1);  % ?????????
   str{sc}=['lmiterm([' nms ' 1 1 P],[eye(' nss ...
         ');zeros(' mms ',' nss ')],state0,''s'');'];
   eval(str{sc});
   sc=sc+1;      
end
for j=1:E.nsimple,
   js=num2str(j);
   nh=0;           % C-position counter for the right (C) terms 
   ng=0;           % C-position counter for the left (B) terms
   for i=1:size(E.X{E.simples(j)},2),
      is=num2str(i);
      vrs=num2str(E.X{E.simples(j)}(1,i));
      eval(['L_' js '_' is ...
            '=E.B{E.simples(j)}(ng+1:ng+E.X{E.simples(j)}(2,i),:)'';'])
      eval(['R_' js '_' is ...
            '=E.gqfm(j)*E.C{E.simples(j)}(nh+1:nh+' ...
            'E.X{E.simples(j)}(3,i),:);']) 
      str{sc}=['lmiterm([' nms ' 1 1 ' vrs '],L_' js '_' is ...
            ',R_' js '_' is ',''s'');'];
      eval(str{sc});
      sc=sc+1;           
      nh=nh+E.X{E.simples(j)}(3,i);
      ng=ng+E.X{E.simples(j)}(2,i);
   end
end
% introduce the gain estimation terms
str{sc}='\n%% L2 gain estimation terms ...';
sc=sc+1;
c_z=E.C{double(z)};
c_f=E.C{double(f)};
str{sc}=['lmiterm([' nms ' 1 1 0],c_z''*c_z)'];
eval(str{sc});
sc=sc+1;
str{sc}='xg=lmivar(1,[1 1]);';
eval(str{sc});
sc=sc+1;
str{sc}=['lmiterm([-' nms ' 1 1 xg],c_f''*c_f,1);'];
eval(str{sc});
sc=sc+1;
 
% solving the system of LMI's 

str{sc}='\n%% Solving the system of LMIs ...';
sc=sc+1;
str{sc}='lmi=getlmis;';
eval(str{sc});
sc=sc+1;
str{sc}='ndec=decnbr(lmi);';
eval(str{sc});
sc=sc+1;
str{sc}='c=zeros(ndec,1);';
eval(str{sc});
sc=sc+1;
str{sc}='nc=decinfo(lmi,xg);';
eval(str{sc});
sc=sc+1;
str{sc}='c(nc)=1;';
eval(str{sc});
sc=sc+1;
options=E.options;
str{sc}='[gain,xopt]=mincx(lmi,c,options);';
sc=sc+1;
str{sc}='gain=sqrt(gain);';

if E.options(5)>100,
   save lmi_exe state* struct* options l_* r_* L_* R_* c_f c_z
   fid=fopen('lmi_exe.m','wt');
   for i=1:sc,
      fprintf(fid,[str{i} '\n']);
   end
   fclose(fid);
end

if vrb,
   disp(['  Solving with ' num2str(ndec) ' decision variables ...'])
end

eval(str{sc-1})
eval(str{sc})

if ~isempty(xopt),
   global ABST
   ABST.xopt=xopt;
   ABST.E=E;
end
