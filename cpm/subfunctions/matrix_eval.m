%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     %%
%%  matrix evaluation  %%
%%                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%

function [Xe]=matrix_eval(A,X,B,coeff)
% function [Xe]=matrix_eval(A,X,B,coeff)
%
% evaluate matric function F(x)=A*X*B
%
% written by cykao@mit.edu  Oct. 14 1998

mX=size(X,1);
nX=size(X,2);
lc=length(coeff);

Xv=zeros(mX,nX);

for flcnt1=1:mX      % index of row of matrix X
    for flcnt2=1:nX  % index of column of matrix X
        nth_var=X(flcnt1,flcnt2);
        if nth_var>lc,
           error('Bad index of independent variables')
        else
           if nth_var>0,
              Xv(flcnt1,flcnt2)=coeff(nth_var);
           elseif nth_var<0
              Xv(flcnt1,flcnt2)=-coeff(-nth_var);
           end
        end
     end
end
Xe=A*Xv*B;
