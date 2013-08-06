function p=abst(a,s)
% function p=abst(a,s)
%
% with no arguments, produces an empty "abst" element
%
% with a single element,
% tries to convert it to the "abst" class type,
%
% with two inputs, s=0,
% creates an abst variable referencing the log entry #a
%
% Written by ameg@mit.edu,  last modified October 13, 1997
% Last modified by cmj on 2013/4/18

superiorto('lti')
superiorto('tf')
superiorto('ss')
superiorto('zpk')
superiorto('double')

global ABST

if isempty(ABST),
    disp_str(12)
end

if nargin==0,  % default definition
    p.n=1;
    p=class(p,'abst');
    return
end

if (nargin==2)&&(s==0),
    if ABST.nlog<a,
        disp_str(35)
    else
        p.n=a;
        p=class(p,'abst');
        return
    end
end

if isempty(a),
    p.n=1;
    p=class(p,'abst');
    return
end

intype=class(a);        % type of a

if isa(a,'abst'),
    p=a;
    return
end

if nargin==1,                % conversion
    for k=1:ABST.next,
        if isa(a,ABST.ext{k,1}),
            ABST.ndat=ABST.ndat+1;            % count the data
            ABST.dat{ABST.ndat}=a;            % store the data
            x=zeros(1,ABST.mlog);             % prepare new log entry
            x(1)=ABST.ext{k,2};               % interior class
            x(2)=size(a,1);       % vertical size
            x(3)=size(a,2);       % horizontal size
            x(4)=-1;              % -1 means "conversion"
            x(5)=k;               % converted from external type ...
            x(6)=ABST.ndat;       % data entry address
            z=abst_alloc(x);
            p.n=z;
            p=class(p,'abst');
            return
        end
    end
    disp_str(36,ABST.name,intype)
end

disp_str(37)
