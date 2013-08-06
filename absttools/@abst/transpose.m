function c=transpose(a)
% function c=transpose(a)
%
% transpose function for the "abst" type
%
%
% Written by ameg@mit.edu,  last modified October 13, 1997
% Last modified by cmj on 2013/4/19

global ABST
% is operation supported ?
if ~isfield(ABST,'transpose'),
    disp_str(14,'transpose')
end

a=abst(a);                % convert to "abst", if necessary
na=double(a);             % reference to a in ABST.log
ca=ABST.log(na,1);        % interior class of a
[va,ha]=size(a);

ocls=ABST.transpose(ca);    % check if "transpose" is allowed
if ocls==0,
    disp_str(40,ABST.cls{ca},'.''')
else
    z=abst_alloc([ocls ha va num_op('transpose') 0 na]);
    c=abst(z,0);
end
