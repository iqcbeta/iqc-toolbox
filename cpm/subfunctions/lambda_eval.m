%%%%%%%%%%%%%%%%%%%%%%%%
%%                    %%
%%  plug in \lambda   %% 
%%                    %%
%%%%%%%%%%%%%%%%%%%%%%%%

function MX=lambda_eval(lambda,nth_lmi)

global L


nterm=length(L.lmi{nth_lmi}.X);
MX=zeros(L.lmi{nth_lmi}.size);
for flcnt2=1:nterm
    G = L.lmi{nth_lmi}.G{flcnt2};
    X = L.lmi{nth_lmi}.X{flcnt2};
    H = L.lmi{nth_lmi}.H{flcnt2};
    if ~isempty(X)
       MX= MX + matrix_eval(G,X,H,lambda);  % non-constant term in LMI
    else
       MX= MX + L.lmi{nth_lmi}.C;           % constant term in LMI 
    end
end
MX=0.5*(MX+MX');   %%% all LMIs "M<0" are interpreted as " M+M'<0 " %%%



