function c=le(a,b)
% function c=le(a,b)
% 
% le function for the "abst" type
% 
%
% Written by ameg@mit.edu,  last modified October 13, 1997
global ABST

% is operation supported ?
if ~isfield(ABST,'le'),
 error('Operation "le" not supported')
end

a=abst(a);                % convert to "abst", if necessary
na=double(a);             % reference to a in ABST.log
ca=ABST.log(na,1);        % interior class of a
[va,ha]=size(a);
b=abst(b);
nb=double(b);
cb=ABST.log(nb,1);
[vb,hb]=size(b);

ocls=ABST.le(ca,cb);    % check if "le" is allowed
if ocls==0,
   error([ABST.cls{ca} '<=' ABST.cls{cb} '  not allowed'])
elseif (va==vb)&(ha==hb),
   z=abst_alloc([ocls va hb num_op('le') 0 na nb]);
   c=abst(z,0);
elseif (va==1)&(ha==1),
   z=abst_alloc([ocls vb hb num_op('le') 0 na nb]);
   c=abst(z,0);
elseif (vb==1)&(hb==1),
   z=abst_alloc([ocls va ha num_op('le') 0 na nb]); 
   c=abst(z,0);
else
   error('argument dimensions are not compatible')
end

   

