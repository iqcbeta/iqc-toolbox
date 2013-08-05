function [gain,xopt]=iqc_gain_accpm(varargin)
% 
% Analytic Center Cutting Plane Method for IQC problems 
%
% function [gain,xopt]=lmi_gain_accpm(f,y,options)
%
% gives best estimate of the L2 gain f->y in the system of IQC's
% defined using the "abst" environment "iqc". f and y can be of type "inp" or "sgn".
%
% The 'options' parameter is optional. This input must be a structure which
% has the following fields: 
%
%   options.max_iteration : Program is forced to stop when number of iteration > this number
%   options.boxlen        : Each decision variable is assumed to be in [-boxlen, +boxlen]  
%   options.accind1       : Program will stop when (obj_current - obj_lowerbound)/ obj_current < accind1, if obj_current > 0.001
%   options.accind2       : Program will stop when (obj_current - obj_lowerbound) < accind2, if obj_current <= 0.001
%   options.accind3       : Program will stop when the volumn of the outbound ellispoid is less than accind3
%
% The default values of
%
%   options.max_iteration = 100000
%   options.boxlen        = 1e7
%   options.accind1       = 1e-4
%   options.accind2       = 1e-6
%   options.accind3       = 1e-12
%
% If the problem has a solution, the iqc_value.m can be used to
% store the optimal multipliers in global variable ABSTSOLUTION
% (so that they can be retrieved using "value.m") 
%
% Last modified by cykao@ee.unimelb.edu.au on June 06 2006 


if nargin == 3,
   if isstruct(varargin{3})~=1, 
      error('The ''options'' must be a structure') 
   end
   f = varargin{1};
   y = varargin{2};
   options = varargin{3};
elseif nargin<3,
   if nargin<2, 
      error('Two input signals are required') 
   end
   f = varargin{1}; 
   y = varargin{2};
   options = [];
else
   error('Too many inputs')
end

% --------------------------------- Parameters ----------------------------------------------------------
if isempty(options) 
   max_iteration = 100000;   
   boxlen        = 1e7;      
   accuracy_ind1 = 1e-4;     
   accuracy_ind2 = 1e-6;
   accuracy_ind3 = 1e-12;
   fset          = {};
else   
   if isfield(options,'max_iteration')
      max_iteration = options.max_interation;
   else
      max_iteration = 100000;
   end
   
   if isfield(options,'boxlen')
      boxlen  = options.boxlen ;
   else
      boxlen  = 1e7;
   end

   if isfield(options,'accind1')
      accuracy_ind1 = options.accind1;
   else
      accuracy_ind1 = 1e-4;
   end
      
   if isfield(options,'accind2')
      accuracy_ind2 = options.accind2;
   else
      accuracy_ind2 = 1e-6;
   end

   if isfield(options,'accind3')
      accuracy_ind3 = options.accind3;
   else
      accuracy_ind3 = 1e-12;
   end
   
   if isfield(options,'fset')
      fset = options.fset;
   else
      fset = {};
   end

end
positive_ind  = 1e-5;
% -------------------------------------------------------------------------------------------------------

global L LP_A LP_B
global ndecv_e ndecv_o
L=iqc_pre_optlp(f,y,1);

if isempty(L.obj)
   feap = 1;
else
   feap = 0;
end

err_msg       = 0;
ndecv_o       = L.ndecvar;
no_feasible   = 1;
reset         = 0;
obj_now       = [];

if isempty(fset)
   A_box  = [eye(ndecv_o); -eye(ndecv_o)];
   B_box  = [boxlen*ones(ndecv_o,1); boxlen*ones(ndecv_o,1)];
else
   if ~isfield(fset,'center') | ~isfield(fset,'range') 
      error('Error in option setting: the initial feasible set ''fset'' needs to have two fields specified: ''center'' and ''range''')
   else
      A_box  = [eye(ndecv_o); -eye(ndecv_o)];
      B_box  = [fset.range + fset.center; fset.range - fset.center];
   end
end

% ------ set options ------- %
%  eig_options.disp    = 0;  % 
% -------------------------- %



% +--------------------------------------------+
% |                                            |
% |  find out the non-frequence depedent LMIs  |
% |  and frequence depedent LMIs               |
% |                                            |
% +--------------------------------------------+
fdeplmi  = [];
nfdeplmi = [];
phc_A    = [];
phc_B    = [];
for flcnt1=1:length(L.lmi)
    if isempty(L.lmi{flcnt1}.A)
       if L.lmi{flcnt1}.size==1,  % collect original polyhedral constraints 
          [coeff,constant]=hyperplane(L,flcnt1,1,ndecv_o);
          norm_cnst = norm(coeff,2); 
          coeff     = real(coeff./norm_cnst);       % normalize the constraint
          constant  = real(-constant/norm_cnst);
          phc_A=[phc_A;coeff   ];
          phc_B=[phc_B;constant];
          L.lmi{flcnt1}.done=1;
       else
          nfdeplmi=[nfdeplmi;flcnt1];
       end
    else
       fdeplmi=[fdeplmi;flcnt1];
    end
end
LP_A = [[zeros(1,ndecv_o-1) -1]; phc_A];
LP_B = [0; phc_B];


% +------------------------+ 
% |                        |
% |       Main Loop        |    
% |                        |
% +------------------------+

done=0;          % flag : set to 1 when optimal solution found
counter=1;       % counter for how many loops the program has been running
counterp=1;
feasible=0;      % feasibility flag; feasible=1 if feasible point found.

disp(' ')
disp('Now, solving LMIs ......')
disp(' ')
%%             1         2         3         4         5         6         7         8
%% ## 1234567890123456789012345678901234567890123456789012345678901234567890123456789012345  
disp(' Iteration    |   Best obj. so far    |   Lower bound of opt. obj.    |     Gap  ')
disp(' -------------+-----------------------+-------------------------------+---------------- ')
disp(' ')

LB = [];

while ~done & counter<=max_iteration,
 
 %% ======================= %%
 %%                         %%
 %%  Generate a test point  %%
 %%                         %%
 %% ======================= %%
 
 go_home = 0;
 
 if counter == 1

    Atot = [A_box; LP_A];
    Btot = [B_box; LP_B];
    [xc, opt, Del2F, Vol1, eta] = cmstep(2, Atot, Btot);
    Vol0    = Vol1;
    DD2F1   = sqrt(det(Del2F));
    DD2F0   = DD2F1;
    nc_curr = size(Atot,1);
     
 else
 
    Atot = [A_box; LP_A];
    Btot = [B_box; LP_B];
    [xc, opt, Del2F, Vol1, eta] = cmstep(2, Atot, Btot, xfeas);
    DD2F1   = sqrt(det(Del2F));
    nc_curr = size(Atot,1);
    
 end

 %% =================== %%
 %%                     %%
 %%   constraint drop   %%
 %%                     %%
 %% =================== %%

 if eta < 0.5         
    r = size(Atot,1)/(1-eta) + eta/(1-eta); 
    
    index_collect_box = [];
    for flcnt = 1:size(A_box,1)
        dist1 = abs(B_box(flcnt) - A_box(flcnt,:)*xc);
        dist2 = sqrt(A_box(flcnt,:)*Del2F*A_box(flcnt,:)')*r;
        if dist1 <= dist2 
           index_collect_box(flcnt) = flcnt;
        end
    end
    index_collect_box = find(index_collect_box > 0); % -- indices of the constraints not to be removed
    
    index_collect_LP = [];
    for flcnt = 1:size(LP_A,1)
        dist1 = abs(LP_B(flcnt) - LP_A(flcnt,:)*xc);
        dist2 = sqrt(LP_A(flcnt,:)*Del2F*LP_A(flcnt,:)')*r;
        if dist1 <= dist2 
           index_collect_LP(flcnt) = flcnt;
        end
    end
    index_collect_LP = find(index_collect_LP > 0); % -- indices of the constraints not to be removed
    
    A_box = A_box(index_collect_box,:);
    B_box = B_box(index_collect_box,:);
    LP_A  = LP_A(index_collect_LP,:);
    LP_B  = LP_B(index_collect_LP,:);
 end
 
 
 %% =================== %%
 %%                     %%
 %%   Revoking Oracle   %%
 %%                     %%
 %% =================== %%


 % ---- 1. check the non-freq depedent LMIs ----

    for flcnt1=1:length(nfdeplmi)
        nth_lmi = nfdeplmi(flcnt1);
        MX = lambda_eval(xc,nth_lmi);
        [cholR,cholp]=chol(-MX);
        if cholp>0,           
           [V,D] = eig(MX);
           [D IND] = max(diag(D));
           V = V(:,IND);               
           if ~(D> -positive_ind)
              if err_msg == 1
                 disp(['Warning (1): A supposed-to-be non-negative definite matrix has maximum eigenvalue ',sft(num2str(D),12)])
              end
           else
              coeff=hyperplane(L,nth_lmi,V,ndecv_o);
              norm_cnst = norm(coeff,2); 
              coeff     = real(coeff./norm_cnst); % normalize the constraint
              constant  = coeff*xc;
              LP_A      = [LP_A; coeff];          % incorperate the constraint with previous constraints
              LP_B      = [LP_B; constant];
              go_home   = 1;
              break
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
               [V,D] = eig(R);
               [D IND] = max(diag(D));
               V = V(:,IND);               
               if ~(D> -positive_ind)
                   if err_msg == 1
                      disp(['Warning (2): A supposed-to-be non-negative definite matrix has maximum eigenvalue ',sft(num2str(D),12)])
                   end
               else
                   vp        = [zeros(nA,1);V];
                   coeff     = hyperplane(L,nth_lmi,vp,ndecv_o);
                   norm_cnst = norm(coeff,2); 
                   coeff     = real(coeff./norm_cnst); % normalize the constraint
                   constant  = coeff*xc;
                   LP_A      = [LP_A; coeff];          % incorperate the constraint with previous constraints
                   LP_B      = [LP_B; constant];
                   go_home   = 1;
                   break
               end

            else 

            % ---------------------------------------------------------------- %
            %   Step 2: <<R is negative definite>> check if the Hamiltonian    %
            %   matrix has eigenvalues on the imaginary axis                   %
            % ---------------------------------------------------------------- %

       	       omega = IS_Hsys_OK(A,B,Q,F,R);
               if ~isempty(omega)
		          found = 0;
		          ctn   = 1;
		          while (~found & ctn <= length(omega))
                         GG1=[inv(j*omega(ctn)*eye(size(A,1))-A)*B;eye(size(B,2))];
                         GG=GG1'*[Q,F;F',R]*GG1;
                         [VGG,DGG] = eig(GG);
                         [DGG IND] = max(real(diag(DGG)));
		                 if ~(DGG> -positive_ind)
			                  DGGS(ctn)= DGG;
		                      VECS{ctn}= GG1*VGG(:,IND);
			             else
			                  VGG = VGG(:,IND);
		                      found = 1;
		                 end
			             ctn = ctn + 1;
	               end

                   if ~found
                      if err_msg == 1
                         disp(['Warning (3): A supposed-to-be non-negative definite matrix has maximum eigenvalue ',sft(num2str(max(DGGS)),12)])
                      end
                   else
                      vec=GG1*VGG;
                      coeff=hyperplane(L,nth_lmi,vec,ndecv_o);
                      norm_cnst = norm(coeff,2); 
                      coeff     = real(coeff./norm_cnst); % normalize the constraint
                      constant  = coeff*xc;
                      LP_A      = [LP_A; coeff];          % incorperate the constraint with previous constraints
                      LP_B      = [LP_B; constant];
                      go_home   = 1;
                      break
                   end
               end
            end
        end

    end


 %% ======================== %%
 %%                          %%
 %%   Oracle say 'yes' ...   %%
 %%                          %%
 %% ======================== %%

 if ~go_home
    if feap 
       no_feasible = 0;
       done = 1;
       xopt = xc;
    else 
       if no_feasible
          no_feasible = 0;
          obj_now = L.obj*xc;
          LP_A  = [ L.obj; LP_A]; 
          LP_B  = [ obj_now; LP_B];
          coeff = L.obj;
          xopt  = xc;
       else
          obj_now = L.obj*xc;
          LP_B(1) = obj_now;
          coeff   = L.obj;
          xopt    = xc; 
       end
    end
 end 


 %% ================== %%
 %%                    %%
 %%  Stoping criteria  %%
 %%                    %%
 %% ================== %%

 if ~feap
    if isempty(obj_now) 
       str1 = sft('------',12); 
       str2 = sft('------',12);
       str3 = sft('------',12);
    else 
       str1 = sft(num2str(obj_now),12);
       if eta < 0.5
          LB = L.obj*xc - (nc_curr/(1-eta)+eta/(1-eta))*sqrt(L.obj*Del2F*L.obj');
          gap  = obj_now - LB;
       end
       
       if ~isempty(LB) 
           str2 = sft(num2str(LB),12);     
           str3 = sft(num2str(gap),12);
       else
           str2 = sft('------',12);
           str3 = sft('------',12);               
       end
    end
 end

 if mod(counter,10)==1
     if ~feap
        str4=[blanks(2),sft(num2str(counterp),5),blanks(12),str1,blanks(15),str2,blanks(17),str3];
        disp(str4)
     else
        disp([blanks(2),sft(num2str(counterp),5),blanks(12), ' ......................... ' ] ) 
     end
     Vol0 = Vol1;
     DD2F0 = DD2F1;
     counterp = counterp + 1;
 end

 if ~feap & ~isempty(obj_now)
     if (abs(obj_now)>1e-3 & abs(gap/obj_now) < accuracy_ind1) | (abs(obj_now)<=1e-3 & gap < accuracy_ind2)
        str4=[blanks(2),sft(num2str(counterp),5),blanks(12),str1,blanks(15),str2,blanks(17),str3]; 
        disp(str4)
        done = 1;
     end
 else
     if  Vol1 < accuracy_ind3
         done = 1;
     end
 end

 %% ============================ %%
 %%                              %%
 %%  computing a feasible point  %% 
 %%  for next loop               %%
 %%                              %%
 %% ============================ %%

    beta  = sqrt(2)/2;
    aD    = Del2F*coeff';
    r     = sqrt(coeff*aD);
    dx    = -(beta/r)*aD;
    xfeas = xc+dx;
    s     = [B_box; LP_B] - [A_box; LP_A]*xfeas;
    if any(s<=0)
       disp(' The program is forced to stop --- matrix very ill-conditioned ')
       done = 1;
    end
    counter = counter+1;  % increase counter

end %%% <--- end of Main Loop!! 

if no_feasible
   gain=[];
   xopt=[];
   disp(' ')
   disp(['No feasible point found after ',num2str(counter-1),' iteration'])
else
   if feap
      disp(' ')
      disp(' ...... feasible point found ')
   else
      gain = sqrt(L.obj*xopt);
   end
end

global ABST
ABST.xopt=xopt;












