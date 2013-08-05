function c=display(a)
% function c=display(a)
% 
% command window output function for the "abst" type
% 
%
% Written by ameg@mit.edu,  last modified October 13, 1997
global ABST

k=double(a);

if isempty(k),
  disp([inputname(1) ' = []'])
  return
end

hdr='#     class     vsize hsize operation   comment    arg1 arg2 arg3 ' ;
disp('===================================================================')
disp(hdr)
disp('-------------------------------------------------------------------')
ss=[sft(k,6) sft(ABST.cls{ABST.log(k,1)},10)];
  ss=[ss sft(ABST.log(k,2),6) sft(ABST.log(k,3),6)];
  ss=[ss sft(num_op(ABST.log(k,4)),12)];
  if ABST.log(k,4)==-1,        % conversion
     ss=[ss sft(ABST.ext{ABST.log(k,5),1},11)];
  elseif ABST.log(k,4)==0,     % definition
     ss=[ss sft(ABST.def{ABST.log(k,5)},11)];
  else
   ss=[ss sft('  -  ',11)];
  end
  ss=[ss sft(ABST.log(k,6),5) sft(ABST.log(k,7),5) sft(ABST.log(k,8),5)];
    disp(ss)
disp('===================================================================')