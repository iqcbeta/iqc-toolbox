function c=ctranspose(a)
% function c=ctranspose(a)
% 
% ctranspose function for the "abst" type
% 
%
% Written by ameg@mit.edu,  last modified October 13, 1997
global ABST

% is operation supported ?
if ~isfield(ABST,'ctranspose'),
 error('Operation "ctranspose" not supported')
end

a=abst(a);                % convert to "abst", if necessary
na=double(a);             % reference to a in ABST.log
ca=ABST.log(na,1);        % interior class of a
[va,ha]=size(a);

ocls=ABST.ctranspose(ca);    % check if "ctranspose" is allowed
if ocls==0,
   error([ABST.cls{ca} '''  not allowed'])
else
   z=abst_alloc([ocls ha va num_op('ctranspose') 0 na]);
   c=abst(z,0);
end

   

