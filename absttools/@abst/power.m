function c=power(a,b)
% function c=power(a,b)
% 
% power function for the "abst" type
% 
%
% Written by ameg@mit.edu,  last modified October 13, 1997
global ABST

% is operation supported ?
if ~isfield(ABST,'power'),
 error('Operation "power" not supported')
end

a=abst(a);                % convert to "abst", if necessary
na=double(a);             % reference to a in ABST.log
ca=ABST.log(na,1);        % interior class of a
[va,ha]=size(a);
b=abst(b);
nb=double(b);
cb=ABST.log(nb,1);
[vb,hb]=size(b);

ocls=ABST.power(ca,cb);    % check if "power" is allowed
if ocls==0,
   error([ABST.cls{ca} '.^' ABST.cls{cb} '  not allowed'])
elseif (va==vb)&(ha==hb),
   z=abst_alloc([ocls va hb num_op('power') 0 na nb]);
   c=abst(z,0);
elseif (va==1)&(ha==1),
   z=abst_alloc([ocls vb hb num_op('power') 0 na nb]);
   c=abst(z,0);
elseif (vb==1)&(hb==1),
   z=abst_alloc([ocls va ha num_op('power') 0 na nb]); 
   c=abst(z,0);
else
   error('argument dimensions are not compatible')
end

   

