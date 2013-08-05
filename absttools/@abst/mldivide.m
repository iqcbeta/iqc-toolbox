function c=mldivide(a,b)
% function c=mldivide(a,b)
% 
% mldivide function for the "abst" type
% 
%
% Written by ameg@mit.edu,  last modified October 13, 1997
global ABST

% is operation supported ?
if ~isfield(ABST,'mldivide'),
 error('Operation "mldivide" not supported')
end

a=abst(a);                % convert to "abst", if necessary
na=double(a);             % reference to a in ABST.log
ca=ABST.log(na,1);        % interior class of a
[va,ha]=size(a);
b=abst(b);
nb=double(b);
cb=ABST.log(nb,1);
[vb,hb]=size(b);

ocls=ABST.mldivide(ca,cb);    % check if "mldivide" is allowed
if ocls==0,
   error([ABST.cls{ca} '\' ABST.cls{cb} '  not allowed'])
elseif (ha==vb)&(va==ha),
   z=abst_alloc([ocls va hb num_op('mldivide') 0 na nb]);
   c=abst(z,0);
elseif (va==1)&(ha==1),
   z=abst_alloc([ocls vb hb num_op('mldivide') 0 na nb]); 
   c=abst(z,0);
elseif (vb==1)&(hb==1)&(va==ha),
   z=abst_alloc([ocls va va num_op('mldivide') 0 na nb]); 
   c=abst(z,0);   
else
   error('argument dimensions are not compatible')
end

   

