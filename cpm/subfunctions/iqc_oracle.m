 function [LP_A_add,LP_B_add,go_home] = iqc_oracle(xc,L,nfdeplmi,fdeplmi,positive_ind)
  
 ndecv_o = L.ndecvar;
 go_home = 0;
 LP_A_add = [];
 LP_B_add = [];
 
 % ---- 1. check the non-freq depedent LMIs ----

    for flcnt1=1:length(nfdeplmi)        
        nth_lmi = nfdeplmi(flcnt1);
        MX = lambda_eval(xc,nth_lmi);
        [cholR,cholp]=chol(-MX);
        if cholp>0,           
           [V,D]=eig(MX);
           D = real(diag(D));
           pos = find(D> -positive_ind);
           if isempty(pos)
              if err_msg == 1
                 disp('something is wrong ... 1')
              end 
           else
              for flcnt2 = 1:length(pos)
                  eV = V(:,pos(flcnt2));
                  [coeff,constant]=hyperplane(L,nth_lmi,eV,ndecv_o);
                  norm_cnst = norm(coeff,2);
                  coeff     = real(coeff./norm_cnst);
                  LP_A_add  = [LP_A_add; coeff]; 
                  LP_B_add  = [LP_B_add; coeff*xc];
              end
              go_home = 1;
           end
        end
    end
 
 % ---- 2. check the freq depedent LMIs ----  

 
    if ~go_home   % i.e., the non-freq LMIs passed
        
       for flcnt1=1:length(fdeplmi)
            nth_lmi = fdeplmi(flcnt1);
            MX = lambda_eval(xc,nth_lmi);
            A  = L.lmi{nth_lmi}.A;
            B  = L.lmi{nth_lmi}.B;
            nA = size(A,1);
            mA = size(A,2);
            mB = size(B,2);
            Q  = MX(1:nA,1:mA);
            F  = MX(1:nA,mA+1:mA+mB);
            R  = MX(nA+1:nA+mB,mA+1:mA+mB);

            % ------------------------------------------- %
            %   Step 1: check if R is negative definite   %
            %   <<  R<0 is a necessary condition !! >>    %
            % ------------------------------------------- %

            [cholR,cholp]=chol(-R);

            if cholp>0,
               [V,D]=eig(R);
               D = real(diag(D));
               pos = find(D> -positive_ind);
               if isempty(pos)
                  if err_msg == 1
                     disp('something is wrong ... 2')
                  end
               else
                  for flcnt2 = 1: length(pos)
                      vp= [zeros(nA,1);V(:,pos(flcnt2))];
                      [coeff,constant]=hyperplane(L,nth_lmi,vp,ndecv_o);
                      norm_cnst = norm(coeff,2); 
                      coeff     = real(coeff./norm_cnst);
                      LP_A_add  = [LP_A_add; coeff   ]; 
                      LP_B_add  = [LP_B_add; coeff*xc];
                  end
                  go_home   = 1;
               end

            else 

            % ---------------------------------------------------------------- %
            %   Step 2: <<R is negative definite>> check if the Hamiltonian    %
            %   matrix has eigenvalues on the imaginary axis                   %
            % ---------------------------------------------------------------- %

	          omega = IS_Hsys_OK(A,B,Q,F,R);
              if ~isempty(omega)
   		         for flcnt2 = 1: length(omega)
                     oga=omega(flcnt2);
                     GG1=[inv(j*oga*eye(size(A,1))-A)*B;eye(size(B,2))];
                     GG=GG1'*[Q,F;F',R]*GG1;
                     [VGG,DGG]=eig(GG);
                     DGG = real(diag(DGG));
                     pos = find(DGG>-positive_ind);
                     if isempty(pos)
                        if err_msg == 1
                           disp('something is wrong ... 3')
                        end
                     else 
                        for flcnt3 = 1: length(pos)
                            vec=GG1*VGG(:,pos(flcnt3));
                           [coeff,constant]=hyperplane(L,nth_lmi,vec,ndecv_o);
                           norm_cnst = norm(coeff,2); 
                           coeff     = real(coeff./norm_cnst);
                           LP_A_add  = [LP_A_add; coeff   ]; 
                           LP_B_add  = [LP_B_add; coeff*xc];
                        end
                     end
                  end
                  go_home = 1;
              end
            end 
       end
    end
