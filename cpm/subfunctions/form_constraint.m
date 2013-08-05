%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     %%
%%   form_constraint   %%
%%                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%

function [coeff]=form_constraint(A,B,X,v,nvar)
% function [coeff]=form_constraint(A,B,X,v,nvar)
% 
% given v, and matrix function A'*X*B, where X is a matrix 
% variable; this program computes the coefficients of the 
% hyperplane H(x)=(v'A')*X*(Bv) 
%
% written by cykao@mit.edu  Oct. 14 1998

mX=size(X,1);
nX=size(X,2);
v_left=A*v;
v_right=B*v;

coeff=zeros(1,nvar);
if (mX==1 & nX==1)
  nth_var=X(1,1);
  if nth_var>0
     coeff(nth_var)=v_left'*v_right;
  elseif nth_var<0
     coeff(-nth_var)=-v_left'*v_right;
  end
else
  for flcnt1=1:mX                    % index of row of matrix X
      for flcnt2=1:nX                % index of column of matrix X
          nth_var=X(flcnt1,flcnt2);  % index of independent variables
          if nth_var~=0,                
             if nth_var>nvar,
                error('Bad index of independent variables!!')
             else
                if nth_var>0
                   coeff(nth_var)=coeff(nth_var)+conj(v_left(flcnt1))*v_right(flcnt2);
                elseif nth_var<0
                   coeff(-nth_var)=coeff(-nth_var)-conj(v_left(flcnt1))*v_right(flcnt2);
                end
             end
          end
      end
  end
end
