function x=vertcat(varargin)
% function x=vertcat(varargin)
%
% vertcat function for the "abst" type
%
%
% Written by cykao@mit.edu and ameg@mit.edu,  last modified October 13, 1997
% Last modified by cmj on 2013/4/19

x=[];
for k=1:size(varargin,2),
    if ~isempty(varargin{k})
        if isempty(x)
            x=varargin{k};
        else
            x=vcat(x,varargin{k});
        end
    end
end

function c=vcat(a,b)

global ABST

% is operation supported ?
if ~isfield(ABST,'vertcat'),
    disp_str(14,'vertcat')
end

a=abst(a);                % convert to "abst", if necessary
na=double(a);             % reference to a in ABST.log
b=abst(b);
nb=double(b);

if isempty(na),
    if isempty(nb)
        c=[];
        return;
    else
        c=b;
    end
elseif isempty(nb),
    c=a;
    return;
else
    ca=ABST.log(na,1);        % interior class of a
    [va,ha]=size(a);
    cb=ABST.log(nb,1);
    [vb,hb]=size(b);
    
    ocls=ABST.vertcat(ca,cb);    % check if "vertcat" is allowed
    if ocls==0,
        disp_str(42,ABST.cls{ca},';',ABST.cls{cb})
    elseif ha==hb,
        z=abst_alloc([ocls va+vb ha num_op('vertcat') 0 na nb]);
        c=abst(z,0);
    else
        disp_str(38)
    end
end
