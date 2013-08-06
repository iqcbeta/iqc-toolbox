function c=subsref(a,s)
% function c=subsref(a,s)
% 
% subsref function for the "abst" type
% 
% for subrefs producing scalars, c=a(i,j), an entry
% [type vsize hsize num_op('subsref') 0 #a i j]
% is added to ABST.log, 
%
% when the result is not a scalar, the re-assignment Nx4 matrix [I J P Q]
% is converted to cst v and stored, refering to c(i,j)=a(p,q)
% the corresponding entry is
% [type vsize hsize num_op('subsref') 0 #a #v 0]
%
% Written by ameg@mit.edu,  last modified October 13, 1997

global ABST

% is operation supported ?
if ~isfield(ABST,'subsref'),
 disp_str(14,'subsref')
end

na=double(a);             % reference to a in ABST.log
ca=ABST.log(na,1);        % interior class of a
[va,ha]=size(a);

ocls=ABST.subsref(ca);    % check if "subsref" is allowed
if ocls==0,
   disp_str(43,ABST.cls{ca},'[s]','')
else
   Mvert=repmat((1:va)',1,ha);
   Mhorz=repmat((1:ha),va,1);
   st='(';
   for k=1:length(s.subs),
      if isstr(s.subs{k}),
         st=[st ':,'];
      else
         st=[st '[' num2str(s.subs{k}) '],'];
      end
   end
   st=[st(1:(length(st)-1)) ');'];
   eval(['Nvert=Mvert' st])
   eval(['Nhorz=Mhorz' st])  
   if length(Nvert)==1,
      z=abst_alloc([ocls 1 1 num_op('subsref') 0 na Nvert Nhorz]);
      c=abst(z,0);
   else
      [I,J]=find(Nvert);
      if size(Nvert,1)==1 
         N=[I' J' Nvert(sub2ind(size(Nvert),I',J'))' ...
         Nhorz(sub2ind(size(Nvert),I',J'))'];
      else
         N=[I J Nvert(sub2ind(size(Nvert),I,J)) ...
         Nhorz(sub2ind(size(Nvert),I,J))];
      end
      v=abst(N);
      nv=double(v);
      z=abst_alloc([ocls size(Nvert) num_op('subsref') 0 na nv]);
      c=abst(z,0);
   end
end

   

