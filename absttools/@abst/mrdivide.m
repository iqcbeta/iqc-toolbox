function c=mrdivide(a,b)
% function c=mrdivide(a,b)
%
% mrdivide function for the "abst" type
%
%
% Written by ameg@mit.edu,  last modified October 13, 1997
% Last modified by cmj on 2013/4/19

global ABST

% is operation supported ?
if ~isfield(ABST,'mrdivide'),
    disp_str(14,'mrdivide')
end

a=abst(a);                % convert to "abst", if necessary
na=double(a);             % reference to a in ABST.log
ca=ABST.log(na,1);        % interior class of a
[va,ha]=size(a);
b=abst(b);
nb=double(b);
cb=ABST.log(nb,1);
[vb,hb]=size(b);

ocls=ABST.mrdivide(ca,cb);    % check if "mrdivide" is allowed
if ocls==0,
    disp_str(39,ABST.cls{ca},'/',ABST.cls{cb})
elseif (ha==vb)&&(vb==hb),
    z=abst_alloc([ocls va hb num_op('mrdivide') 0 na nb]);
    c=abst(z,0);
elseif (vb==1)&&(hb==1),
    z=abst_alloc([ocls va ha num_op('mrdivide') 0 na nb]);
    c=abst(z,0);
elseif (va==1)&&(ha==1)&&(vb==hb),
    z=abst_alloc([ocls vb vb num_op('mrdivide') 0 na nb]);
    c=abst(z,0);
else
    disp_str(38)
end
