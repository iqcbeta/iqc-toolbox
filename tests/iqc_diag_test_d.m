A=[-8.9130,1.5647; 2.5647,-3.1850];
B=eye(2);
C=2*eye(2);
G=ss(A,B,C,zeros(2));
G=c2d(G,1);

abst_init_iqc(1)
f=signal(2);
w=signal(2);
v=G*(f+w);
w==iqc_diag(v); %#ok<*EQEFF>
gain=iqc_gain_tbx(f,v)
disp('the result was gain=3.6437 on may 5, 2013');