A0=[0 -1;1 -0.5];
B0=[-2 0;1 0];
C0=[-100 100];
D0=[0 1];
L0=[1 0];
delta=0.3;

nA=size(A0,1);

abst_init_iqc

w=signal(size(B0,2));

G0=ss(A0,eye(nA),eye(nA),zeros(nA,nA));

p=signal;

x=G0*(B0*w+[delta;0]*p);
q=[0 1]*x;
v=L0*x;
y=C0*x+D0*w;

p==iqc_ltigain(q,0.05); %#ok<*EQEFF>

W=1;

% setlmioptions('yalmip','solver','sdpt3')
setlmioptions('lmilab')

[design_min_gamma,estimation]=...
    solveiqc('estimation','l2gain-eli',{w,v,y,p,q,W});

if ~isempty(estimation)
    z=W*(v-estimation*y);
    setlmioptions('lmilab')
    analysis_min_gamma=solveiqc('analysis','l2gain',{w,z});
    disp(['design_l2gain: ',num2str(design_min_gamma)])
    disp(['analysis_l2gain: ',num2str(analysis_min_gamma)])
end