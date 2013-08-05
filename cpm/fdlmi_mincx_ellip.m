function [cost,xopt]=fdlmi_mincx_ellip(varargin)
% 
% Ellipsoid Method for solving frequency dependent LMI
%
% function [cost,xopt]=fdlmi_mincx_ellip(obj,options)
%
% This program solves the LMI problem:
% 
%       Minimize    c'x   subject to    L(w,x)  <  0
%       where L(w,x) is a matrix function, which can be frequency dependent.
%       If L is frequency dependent, then L(w,x) < 0 should be understood
%       as v'Re(L(w,x))v < 0 for all w, and for all v in R^n (w is the frequency) 
%
%  The objective c'x (obj) has to be defined using the "abst" environment "fdlmi". 
%  If the problem is feasible, a global cell array ABSTSOLUTION, of same size as ABST.log, 
%  can be created using iqc_value.m to contain the optimal values of the decision variables. 
%  (so that they can be retrieved using "value.m")
%
% The 'options' parameter is optional. This input must be a structure which
% has the following fields: 
%
%   options.max_iteration : Program is forced to stop when number of iteration > this number
%   options.Rad           : Feasible set is assumed to be contained in the ball of radius = Rad  
%   options.accind1       : Program will stop when (obj_current - obj_lowerbound)/ obj_current < accind1, if obj_current > 0.001
%   options.accind2       : Program will stop when (obj_current - obj_lowerbound) < accind2, if obj_current <= 0.001
%   options.accind3       : Program will stop when the volumn of the outbound ellispoid is less than accind3
%
% The default values of
%
%   options.max_iteration = 100000
%   options.Rad           = 1e7
%   options.accind1       = 1e-4
%   options.accind2       = 1e-6
%   options.accind3       = 1e-12
%
% Last modified by cykao@ee.unimelb.edu.au on June 06 2006 

if nargin == 2,
   if isstruct(varargin{2})~=1, 
      error('The ''options'' must be a structure') 
   end
   obj = varargin{1};
   options = varargin{2};
elseif nargin == 1,
   obj = varargin{1};
   options = [];
elseif nargin == 0, 
   obj = [];
   options = [];
else
   error('Too many inputs')
end

% --------------------------------- Parameters ----------------------------------------------------------
if isempty(options) 
   max_iteration = 100000;   
   Rad           = 1e7;      
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
   
   if isfield(options,'Rad')
      Rad  = options.Rad;
   else
      Rad  = 1e7;
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

% +----------------+
% | initialization |
% +----------------+
if nargin==0,
   E=fdlmi_extract;
   lmi_lp_pre1(E);
else
   E=fdlmi_extract(obj);
   lmi_lp_pre1(E,obj);
end

global L LP_A LP_B
global ndecv_e ndecv_o

if isempty(L.obj)
   feap = 1;
else
   feap = 0;
end

err_msg       = 0;
ndecv_o       = L.ndecvar;
dim           = ndecv_o;
no_feasible   = 1;
obj_now       = [];
phc_A         = [];
phc_B         = [];

% ------ set options ------- %
%  eig_options.disp    = 0;  %
% -------------------------- %

if dim<2
   error('By some stupid reason, this program currently does not work for the situation when there is only 1 decision variable')
end


% +--------------------------------------------+
% |                                            |
% |  find out the non-frequence depedent LMIs  |
% |  and frequence depedent LMIs               |
% |                                            |
% +--------------------------------------------+
fdeplmi  = [];
nfdeplmi = [];
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
LP_A = phc_A;
LP_B = phc_B;


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


while ~done & counter<=max_iteration,
 
 %% ======================= %%
 %%                         %%
 %%  Generate a test point  %%
 %%                         %%
 %% ======================= %%
 

 if counter == 1
    if isempty(fset)
       ELLIP = (Rad)^2*eye(dim);
       xc    = zeros(dim,1);
    else
       if ~isfield(fset,'center') | ~isfield(fset,'range') 
          error('Error in option setting: the initial feasible set ''fset'' needs to have two fields specified: ''center'' and ''range''')
       else
          ELLIP = (fset.range)^2*eye(dim);
          xc    = fset.center;
       end
    end
    
    if ~isempty(LP_A)
       Ind   = find(LP_B - LP_A*xc >=0);
       if ~isempty(Ind)
           [xc, ELLIP, Vol1] = cmstep(3,xc,ELLIP,LP_A,LP_B,Ind);
       else
           Vol1   = pi^(dim/2)*Rad^(dim)/gamma(dim/2+1);
       end
    else
       Vol1   = pi^(dim/2)*Rad^(dim)/gamma(dim/2+1);
    end
    Vol0 = Vol1;
 else
    [xc, ELLIP, Vol1] = cmstep(3,xc,ELLIP,LP_A,LP_B,Ind); 
 end

 if isempty(xc)
    disp(' The program is forced to stop --- matrix very ill-conditioned ')
    break
 end


 %% =================== %%
 %%                     %%
 %%   Revoking Oracle   %%
 %%                     %%
 %% =================== %%

 go_home = 0;
 Ind      = [];
 LP_ADD_A = [];
 LP_ADD_B = [];

 % ---- 0. check the linear constraints ----
 
 if ~isempty(LP_A)
    Ind = find(LP_A*xc - LP_B >= 0);
    if ~isempty(Ind)    
       go_home = 1;
    end
 end

 % ---- 1. check the non-freq depedent LMIs ----
 
 if ~go_home % check Matrix inequalities only when linear constraints are all satisfied
    
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
                 disp(['Warning (1): A supposed-to-be non-negative definite matrix has maximum eigenvalue ',sft(num2str(max(D)),12)])
	      end
           else
              for flcnt2 = 1:length(pos)
                  eV = V(:,pos(flcnt2));
                  [coeff,constant]=hyperplane(L,nth_lmi,eV,ndecv_o);
                  norm_cnst = norm(coeff,2);
                  coeff     = real(coeff./norm_cnst);
                  LP_ADD_A  = [LP_ADD_A; coeff]; 
                  %LP_ADD_B  = [LP_ADD_B; coeff*xc];
                  LP_ADD_B  = [LP_ADD_B; real(-constant/norm_cnst)];
              end
              go_home = 1;
           end
        end
    end
 
 % ---- 2. check the freq depedent LMIs ----  

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
                     disp(['Warning (2): A supposed-to-be non-negative definite matrix has maximum eigenvalue ',sft(num2str(max(D)),12)])
		  end
               else
                  for flcnt2 = 1: length(pos)
                      vp= [zeros(nA,1);V(:,pos(flcnt2))];
                      [coeff,constant]=hyperplane(L,nth_lmi,vp,ndecv_o);
                      norm_cnst = norm(coeff,2); 
                      coeff     = real(coeff./norm_cnst);
                      LP_ADD_A  = [LP_ADD_A; coeff   ]; 
                      %LP_ADD_B  = [LP_ADD_B; coeff*xc];
                      LP_ADD_B  = [LP_ADD_B; real(-constant/norm_cnst) ];
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
               go_home = 1; 
   	       for flcnt2 = 1: length(omega)
                   oga=omega(flcnt2);
                   GG1=[inv(j*oga*eye(size(A,1))-A)*B;eye(size(B,2))];
                   GG=GG1'*[Q,F;F',R]*GG1;
                   [VGG,DGG]=eig(GG);
                   DGG = real(diag(DGG));
                   pos = find(DGG>-positive_ind);
                   if isempty(pos)
		              if err_msg == 1
                         disp(['Warning (3): A supposed-to-be non-negative definite matrix has maximum eigenvalue ',sft(num2str(max(DGG)),12)])
	                  end
                   else 
                      for flcnt3 = 1: length(pos)
                          vec=GG1*VGG(:,pos(flcnt3));
                          [coeff,constant]=hyperplane(L,nth_lmi,vec,ndecv_o);
                          norm_cnst = norm(coeff,2); 
                          coeff     = real(coeff./norm_cnst);
                          LP_ADD_A  = [LP_ADD_A; coeff   ]; 
                          %LP_ADD_B  = [LP_ADD_B; coeff*xc];
                          LP_ADD_B  = [LP_ADD_B; real(-constant/norm_cnst)];
                      end
                   end    
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
          LP_A    = [ L.obj; LP_A]; 
          LP_B    = [ obj_now; LP_B];
       else
          obj_now = L.obj*xc;
          LP_B(1) = obj_now;          
       end
       Ind  = 1;
       xopt = xc;
    end   
 elseif isempty(Ind)
    Ind  = [size(LP_A,1)+1:1:size(LP_A,1)+size(LP_ADD_A,1)];
    LP_A = [LP_A; LP_ADD_A];
    LP_B = [LP_B; LP_ADD_B];
 end 


 %% ================== %%
 %%                    %%
 %%  Stoping criteria  %%
 %%                    %%
 %% ================== %%

 if ~feap       
    if isempty(obj_now) 
       str1 = sft('-------',12); 
       str2 = sft('-------',12);
       str3 = sft('-------',12);
    else 
       str1 = sft(num2str(obj_now),12); 
       LB   = L.obj*xc - sqrt(L.obj*ELLIP*L.obj');%'
       gap  = obj_now - LB;
       str2 = sft(num2str(LB),12);
       str3 = sft(num2str(gap),12);
    end
 end

  if mod(counter,10) == 1 
     if ~feap
        str4=[blanks(2),sft(num2str(counterp),5),blanks(12),str1,blanks(15),str2,blanks(17),str3];
        disp(str4)
     else
        disp([blanks(2),sft(num2str(counterp),5),blanks(12), ' ......................... ']  ) 
     end
     Vol0 = Vol1;
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

  counter=counter+1;  % increase counter

end %%% <--- end of Main Loop!! 

if no_feasible
   cost=[];
   xopt=[];
   disp(' ')
   disp(['No feasible point found after ',num2str(counter-1),' iteration'])
else
   if feap
      cost = [];
      disp(' ')
      disp(' ...... feasible point found ')
   else
      cost = L.obj*xopt;
   end
end

global ABST
ABST.E=E;
ABST.xopt=xopt;
