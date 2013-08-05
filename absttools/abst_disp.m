function abst_disp(arg)
% function abst_disp(arg)
%
% This program is designed to display the log of abst environment  
% 
% 1. display the whole log:
%    This is the default function, no input argument 'arg'
%
% 2. display specific log entries:
%    'arg' = column or row vector containing the log entries which
%            is about to be displayed.
%
%    Ex: iqc_displog([1:10])  :  display 1st to 10th log entries
%        iqc_displog([1 3 4]) :  display 1st, 3rd, 4th log entries
%
% 3. display log entries of specific class of abst objects:
%    'arg' = name of the abst object.
%
%    Ex: iqc_displog('variable') : display all log entries about
%                                  variables
%
% written by ameg@mit.edu   Oct. 13, 1997
% revised by cykao@mit.edu  Oct. 11, 1998
%            last modified  Nov. 20, 1998

global ABST

cst=1;
var=2;
lin=3;
lmi=4;
inp=5;
sgn=6;
vsg=7;
csg=8;
vcs=9;
qfm=10;
iqc=11;
lnk=12;

if isempty(ABST),
   disp(' ')
   disp('The "abst" environment not defined')
else
   if nargin==0,
      entry_disp=1:ABST.nlog;
   else
      switch class(arg)
        case 'char'
             badname=1;
             for flcnt1=1:ABST.ncls
                 if strcmp(arg,ABST.cls{flcnt1})
                    arg=find(ABST.log(:,1)==flcnt1);
                    badname=0;
                    break;
                 end
             end
             if badname,
                error('Bad class name')
             elseif isempty(arg)
                disp(' ')
                disp('There is no log entry for the specified class')
                disp(' ')
                return 
             end
        case 'double'
             arg=arg;
      end
      if max(abs(arg))>ABST.nlog
         error(['Bad index! Log has only #',num2str(ABST.nlog),' entries.'])
      else
         entry_disp=round(abs(arg));
      end
   end
   disp(' ')
   disp('---------------------------------------------------------------------')
   disp(['                              ' ABST.name ])
   disp('---------------------------------------------------------------------')
   disp(' ')
   %             1         2         3         4         5       
   %    123456789012345678901234567890123456789012345678901234567890
   hdr='#     class     vsize hsize operation   comment    arg1 arg2 arg3' ;
   disp(hdr)
   disp(' ')
   for flcnt=1:length(entry_disp),
       k=entry_disp(flcnt);
       ss=[sft(k,6) sft(ABST.cls{ABST.log(k,1)},10)];
       ss=[ss sft(ABST.log(k,2),6) sft(ABST.log(k,3),6)];
       ss=[ss sft(num_op(ABST.log(k,4)),12)];
       if  ABST.log(k,4)==-1,       % conversion
           ss=[ss sft(ABST.ext{ABST.log(k,5),1},11)];
       elseif ABST.log(k,4)==0,     % definition   
           ss=[ss sft(ABST.def{ABST.log(k,5)},11)];
       else
           ss=[ss sft('  -  ',11)];
       end
       ss=[ss  sft(ABST.log(k,6),5) sft(ABST.log(k,7),5)];
       ss=[ss sft(ABST.log(k,8),8)];
       disp(ss)
   end
end
disp(' ')