function c=or(a,b)
% function c=or(a,b)
%
% or function for the "abst" type
%
%
% Written by ameg@mit.edu,  last modified October 13, 1997
% Last modified by cmj on 2013/4/19

global ABST
% is operation supported ?
if ~isfield(ABST,'or')
    disp_str(14,'or')
end

a=abst(a);                % convert to "abst", if necessary
na=double(a);             % reference to a in ABST.log
ca=ABST.log(na,1);        % interior class of a
[va,ha]=size(a);
b=abst(b);
nb=double(b);
cb=ABST.log(nb,1);
[vb,hb]=size(b);

ocls=ABST.or(ca,cb);    % check if "or" is allowed
if ocls==0,
    disp_str(39,ABST.cls{ca},'|',SBST.cls{cb})
elseif (va==vb)&&(ha==hb),
    z=abst_alloc([ocls va hb num_op('or') 0 na nb]);
    c=abst(z,0);
elseif (va==1)&&(ha==1),
    z=abst_alloc([ocls vb hb num_op('or') 0 na nb]);
    c=abst(z,0);
elseif (vb==1)&&(hb==1),
    z=abst_alloc([ocls va ha num_op('or') 0 na nb]);
    c=abst(z,0);
else
    disp_str(38)
end
