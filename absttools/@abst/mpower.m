function c=mpower(a,b)
% function c=mpower(a,b)
% 
% mpower function for the "abst" type
% 
%
% Written by ameg@mit.edu,  last modified October 13, 1997
global ABST

% is operation supported ?
if ~isfield(ABST,'mpower'),
 error('Operation "mpower" not supported')
end

a=abst(a);                % convert to "abst", if necessary
na=double(a);             % reference to a in ABST.log
ca=ABST.log(na,1);        % interior class of a
[va,ha]=size(a);
b=abst(b);
nb=double(b);
cb=ABST.log(nb,1);
[vb,hb]=size(b);

ocls=ABST.mpower(ca,cb);    % check if "mpower" is allowed
if ocls==0,
   error([ABST.cls{ca} '^' ABST.cls{cb} '  not allowed'])
elseif (ha==va)&(hb==1)&(vb==1),
   z=abst_alloc([ocls va va num_op('mpower') 0 na nb]);
   c=abst(z,0);
elseif (va==1)&(ha==1)&(vb==hb),
   z=abst_alloc([ocls vb vb num_op('mpower') 0 na nb]); 
   c=abst(z,0);
else
   error('argument dimensions are not compatible')
end

   

