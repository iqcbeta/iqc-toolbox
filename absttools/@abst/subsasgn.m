function c=subsasgn(a,s,b)
% function c=subsasgn(a,s,b)
%
% subsasgn function for the "abst" type
%
% for assignments of a scalar b, a(s)=b, entries
% [type vsize hsize subsasgn j #a #b i]
% are added to ABST.log, referring to a(i,j)=b
%
% for assignments of non-scalars, the re-assignment Nx4 matrix v=[I J P Q]
% such that a(i,j)=b(p,q),
% is converted to cst and stored
%
% Written by ameg@mit.edu,  last modified October 13, 1997
% Last modified by cmj on 2013/4/19

global ABST
% is operation supported ?
if ~isfield(ABST,'subsasgn'),
    disp_str(14,'subsasgn')
end

a=abst(a);
na=double(a);             % reference to a in ABST.log
ca=ABST.log(na,1);        % interior class of a
[va,ha]=size(a);

b=abst(b);
nb=double(b);             % reference to a in ABST.log
cb=ABST.log(nb,1);        % interior class of a
[vb,hb]=size(b);


ocls=ABST.subsasgn(ca,cb);  % check if "subsasgn" is allowed
if ocls==0,
    disp_str(43,ABST.cls{ca},'[s]=',ABST.cls{cb})
else
    Mvert=repmat((1:vb)',1,hb);
    Mhorz=repmat((1:hb),vb,1);
    Nvert=zeros(va,ha);
    Nhorz=zeros(va,ha);
    st='(';
    for k=1:length(s.subs),
        if isstr(s.subs{k}),
            st=[st ':,'];
        else
            st=[st '[' num2str(s.subs{k}) '],'];
        end
    end
    st=[st(1:(length(st)-1)) ')'];
    eval(['Nvert' st '=Mvert;']);
    eval(['Nhorz' st '=Mhorz;']);
    [I,J,S]=find(Nvert);
    if length(I)==1,
        z=abst_alloc([ocls size(Nvert) num_op('subsasgn') J na nb I]);
        c=abst(z,0);
    else
        %      N=[I J Nvert(sub2ind(size(Nvert),I,J)) ...
        %        Nhorz(sub2ind(size(Nvert),I,J))];
        if size(Nvert,1)==1
            N=[I' J' Nvert(sub2ind(size(Nvert),I',J'))' ...
                Nhorz(sub2ind(size(Nvert),I',J'))'];
        else
            N=[I J Nvert(sub2ind(size(Nvert),I,J)) ...
                Nhorz(sub2ind(size(Nvert),I,J))];
        end
        v=abst(N);
        nv=double(v);
        z=abst_alloc([ocls size(Nvert) num_op('subsasgn') 0 na nb nv]);
        c=abst(z,0);
    end
end

%   M=(0:(vb-1))'*hb*ones(1,hb)+ones(vb,1)*(1:hb);
%   N=zeros(va,ha);
%   st='N(';
%   for k=1:length(s.subs),
%      if isstr(s.subs{k}),
%         st=[st ':,'];
%      else
%         st=[st '[' num2str(s.subs{k}) '],'];
%      end
%   end
%   st=[st(1:(length(st)-1)) ')=M;'];
%   eval(st);
%   v=abst(N);
%   nv=double(v);
%   z=abst_alloc([ocls va ha num_op('subsasgn') 0 na nb nv]);
%   c=abst(z,0);
%end
