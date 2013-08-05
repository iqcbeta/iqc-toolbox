%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                   %%
%%  compute ployhegral constraints   %%
%%                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coeff,constant]=hyperplane(L,nth_lmi,vp,ndecv_o)
% compute polyhedral constraint
nterm=size(L.lmi{nth_lmi}.X,2);
coeff=zeros(1,ndecv_o);              % compute constraint 
constant=0;
for flcnt3=1:nterm
    G=L.lmi{nth_lmi}.G{flcnt3}';
    X=L.lmi{nth_lmi}.X{flcnt3};
    H=L.lmi{nth_lmi}.H{flcnt3};
    if ~isempty(X)
       coeff=coeff+form_constraint(G,H,X,vp,ndecv_o);   % non-constant term
    else
       constant=constant+vp'*L.lmi{nth_lmi}.C*vp;       % constant term
    end
end
