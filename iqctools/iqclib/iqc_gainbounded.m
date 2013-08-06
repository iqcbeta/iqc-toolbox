function [w,M]=iqc_gainbounded(v,a)
w=signal(size(v,1));
M=symmetric;
M>0;
v'*(a^2*M)*v>w'*M*w;