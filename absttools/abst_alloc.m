function y=abst_alloc(x)
% function abst_alloc(x)
% 
% writes "x" as a new row to ABST.log
% y is the number this row
% pads "x" with additional zeros, if necessary
% allocates more space for ABST.log, if necessary
%
% Written by ameg@mit.edu,  last modified October 13, 1997

global ABST

if ABST.mlog<length(x),
   error('input argument too long')
else
   ABST.nlog=ABST.nlog+1;
   ABST.log(ABST.nlog,1:length(x))=x;
end

if ABST.nlog==size(ABST.log,1),
 ABST.log=[ABST.log;zeros(100,ABST.mlog)];
end              % increase log size, if necessary

y=ABST.nlog;