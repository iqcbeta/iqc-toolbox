function [w]=iqc_tvscalar(v,D)
% function [w]=iqc_tvscalarV(v,D)
% adds iqc-terms describing the block
% out=d(t)*in, -D<=d(t)<=D,
% in,out must be IQC output numbers (must be same size outputs)
%
% created by Prof. A Megretski, ameg@mit.edu
% modified for new IQC toolbox by C. Kao cykao@mit.edu, 08/31/1997
% Last modified by cmj on 2013/5/5

if nargin<1
    disp_str(1)
end

if nargin<2
   D=1;
end

m=size(v,1);    % dimention of the input signal
w=signal(m);    % define the output signal
X=symmetric(m);
X>0;

% define the IQC of the time varying scalar nonlinear block
if m==1,
   v'*D^2*(X*v)-w'*(X*w)>0;
else
   Y=skew(m);
   v'*D^2*(X*v)+w'*(Y'*v)+v'*(Y*w)-w'*(X*w)>0;
end
