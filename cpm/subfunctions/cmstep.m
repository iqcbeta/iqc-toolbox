function [x1,varargout] = cmstep(options,varargin)
%
% function [x1,varargout] = cmstep(options,varargin)
%
% options = 1 :
%
%    [x1, L1, Vol] = cmstep(options,x0,L0,A0,B0)
%    
%       computing the minimum-volume ellipsoid that contain 
%       the intersection of P and E, where
%           P := { x | A0 x - B0 < 0}
%           E := { x | (x-x0))' L0^{-1} (x-x0) < 1 }. 
%
%       ***** NOTICE : P is assumed to be a half plane, i.e., B0 - A0*x 
%             is assumed to be a scalar. If B0 - A0*x is not a scalar, the 
%             program will automatically use the first component to define P.
%       
%       Inputs  : x0  -- nx1 vector, the center of the ellipsoid E.
%                 L0  -- nxn positive definite matrix that defines E.
%                 A0  -- 1xn vector. A0 x - B0 < 0 defines the half plane P.
%                 B0  -- (optional) a constant. If B0 is missing, then
%                        P is understood as { x | A0(x - x0) < 0 }. 
%
%       Outputs : x1  -- nx1 vector, the center of the minimum-volume ellipsoid.
%                 L1  -- nxn positive definite matrix that defines the 
%                        ellipsoid { x | (x-x1))' L1^{-1} (x-x1) < 1 }   
%                 Vol -- The volume of the ellipsoid { x | (x-x1))' L1^{-1} (x-x1) < 1 }
%
% options = 2 :
%
%    [x1, obj, Del2F, Vol, eta, xp, s] = cmstep(options,A0,B0,xfeas) 
%
%        computing an approximation to the  analytical center, xc, of a bounded polytope 
%        P := { x | A0 x - B0 < 0}. 
%
%        The analytical center of P is the minimizer of the logarithm barrier function 
%        F(x) = - sum_i log( B0(i)-A0(i,j)x(j) ) 
%
%        ***** NOTICE : Make sure that P is bounded. The program is guaranteed
%              to crush if otherwise ... 
%
%        Inputs  : A0    -- mxn matrix. 
%                  B0    -- mx1 matrix.
%                  xfeas -- (optional) a feasible point to P = {x | A0 x - B0 < 0}.
%                          
%        Outputs : x1    -- an approximation to xc. del(F)'*(del^2(F))^{-1}*del(F) < 1e-6
%                           where del(F) is the gradient of F at x1, del^2(F) is the 
%                           Hessian of F at x1.
%                  obj   -- F(x1);
%                  Del2F -- The inverse of the Hessian of F at x1.                  
%                  Vol   -- The volume of the ellipsoid E:= {x | (x-x1)'*Del2F^{-1}*(x-x1)<m^2}.
%                           Notice that P is contained in E. 
%                  eta   -- del(F)'*(del^2(F))^{-1}*del(F). 
%                   xp   -- xp > 0, and A0' * xp = 0. <-- "primal solution"
%                    s   -- B0 - A0*xc

switch options

case 1

% ------------------------ Case 1. ---------------------
% Let 
% P = { x | A0 x < B0 }, Ph = { x | A0 x = B0 }
% E = { x | (x-x0))' L0^{-1} (x-x0) < 1 }
% 
% Assuming Ph is a hyperplane, the following codes 
% generate the minimum-volume ellipsoid which contains
% the intersection of P and E
% ------------------------------------------------------''

  if nargin < 4
     error('Not enough input arguments ...')
  end

  x0  = varargin{1};
  L0  = varargin{2};
  A0  = varargin{3}(1,:)';
  dim = length(x0);

  if nargin >= 5
     B0 = varargin{4}(1);
     A0 = [A0; B0]; 
     [x1, L1] = e_cut(x0,L0,A0);
  else
     [x1, L1] = e_cut(x0,L0,A0);
  end

  varargout(1) = {L1};
  %varargout(2) = {pi^(dim/2)*(sqrt(abs(det(L1)))/gamma(dim/2+1))};
  varargout(2) = {pi^(dim/2)*sqrt(abs(det( L1/gamma(dim/2+1)^(2/dim) )))};
case 2

% ------------------------ Case 2. ---------------------
% Let P = { x | A x < b } 
% 
% Assuming P is bounded, the following codes generate an 
% approximate analytical center of the polytope P.
% ------------------------------------------------------

   if nargin < 3 
     error('Not enough input arguments ...')
   end

   A = varargin{1};
   b = varargin{2};
   
   if nargin >= 4 
      xfeas = varargin{3};
      xfeas_avail = 1;
   else
      xfeas_avail = 0;
   end

   if ~xfeas_avail
   
   % ---- find a feasible point ----
      nA   = size(A,1);
      mA   = size(A,2);
      beta = 0.5;
      x    = zeros(mA,1);
      pos  = find(b<=1);

      bmod      = b;
      bmod(pos) = 1;

      if ~isempty(pos) 
         done = 0;
      else 
         done = 1;
      end
      
      while ~done
        [xa,count,eta,MM] = ana(A,bmod,x);
        brep = max(b(pos), beta*bmod(pos) + (1-beta)*(A(pos,:)*xa));
        bmod(pos) = brep;
        x = xa;
        if all(bmod == b)
           done = 1;
        end
      end
      [xc,count,eta,Del2F,AS,sinv] = ana(A,b,x);
   else
      [xc,count,eta,Del2F,AS,sinv] = ana(A,b,xfeas);
   end
   
   x_primal = sinv - sinv.*sum(AS'*Del2F*AS,2);
   
   mA  = size(A,1);
   nA  = size(A,2);
   E   = b-A*xc;
   x1  = xc;
   dim = size(Del2F,1);
   varargout(1) = {-sum(log(E))};
   varargout(2) = {Del2F};
   varargout(3) = { pi^(nA/2)*(mA^nA)* sqrt(abs(det( Del2F/gamma(nA/2+1)^(2/dim) ))) };
   varargout(4) = {sqrt(eta)};
   varargout(5) = {x_primal};
   varargout(6) = {b-A*xc};

case 3

% ------------------------ Case 3. ---------------------
%
%
%
%
% ------------------------------------------------------

  if nargin < 6
     error('Not enough input arguments ...')
  end

  x0  = varargin{1};
  L0  = varargin{2};
  A0  = varargin{3};
  B0  = varargin{4};
  Ind = varargin{5};
  dim = length(x0);
  Volorg = pi^(dim/2)*(sqrt(abs(det(L0)))/gamma(dim/2+1));
  x1 = [];
  L1 = [];
  
  for flcnt1 = 1: length(Ind)
      g   = A0(Ind(flcnt1),:);
      gb  = B0(Ind(flcnt1));
      [x, L, badflag] = e_cut(x0,L0,[g';gb]);
      if badflag == 0
         [R,p] = chol(L);
         if p == 0
            x1 = x;
            L1 = L;
         end
      end
  end
 
  if ~isempty(x1)
     x0   = x1;
     L0   = L1;
     done = 0;
  else
     x0   = [];
     L0   = [];
     done = 1; 
  end

  while ~done
    cor = find(A0*x0-B0 >= 0);  
    x1 = [];
    L1 = [];
    if ~isempty(cor)
        for flcnt1 = 1: length(cor)
            g   = A0(cor(flcnt1),:);
            gb  = B0(cor(flcnt1));
            [x, L, badflag] = e_cut(x0,L0,[g';gb]);
            if badflag == 0
               [R,p] = chol(L);
               if p == 0
                  x1 = x;
                  L1 = L;
               end
            end
        end
 
        if ~isempty(x1)
            x0   = x1;
            L0   = L1;
        else
            done = 1; 
        end
    else    
        done = 1;
    end
  end

  x1 = x0;
  varargout(1) = { L0 };
  if ~isempty(L0)
      varargout(2) = { pi^(dim/2)*sqrt( det(L0/gamma(dim/2+1)^(2/dim)) ) };
  else
      varargout(2) = { 0 };
  end

end

function [xa,count,eta,Del2F,AS,sinv] = ana(A,b,x) 

s = b-A*x;
if any(s<=0)
   error('xfeas is NOT fesaible') 
end

done = 0;
count = 1;
eta = [];
   
while ~done
   sinv   = 1./s;
   AS     = zeros(size(A,2),length(sinv));
   for flcnt = 1:length(sinv)
       AS(:,flcnt) = A(flcnt,:)'*sinv(flcnt);
   end

   MM    = AS*AS';

% ----------- trial --------------
%   [VMM, DMM] = eig(MM);
%   DMM   = diag(DMM);
%   invMM = VMM*diag(1./DMM)*VMM';
% ---------------------------------
%   ND  = -invMM*AS;
   ND  = -MM\AS;
   dx  = sum(ND,2);
   ds  = -sum(A*ND,2);
   eta = sinv'*ds;

   if abs(eta)<1e-3 
      done = 1;
   else
      if eta<1
         s_next = s + ds;
         x_next = x + dx;
      else
         s_next = s + ds./(1+eta); 
         x_next = x + dx./(1+eta);
      end
    end
   
   if ~done
      if any((b - A*x_next)<0) | any(s_next<0)
         done = 1;
      else
         s = s_next;
         x = x_next;
      end
   end
   count = count + 1;
end
xa    = x;
sinv  = 1./s;
Del2F = MM\eye(size(MM));



function [x1, L1, bad] = e_cut(x0,L0,A0)

dim = length(x0);
if dim ~= length(A0)

     B0      = A0(length(A0));
     A0      = A0(1:length(A0)-1,1);
     A0hat   = L0*A0;
     A0_N    = A0'*A0hat;
     A0_N2   = sqrt(A0'*A0hat);
     AdM     = A0hat*A0hat';
 
     cM      = A0'*x0 + A0_N2/dim;
     cm      = A0'*x0 - A0_N2;

     if B0 >= cM    % The minimum volume ellipsoid is the original one
        x1  = x0;
        L1  = L0;
        bad = 1;
     elseif B0 < cm % The intersection of P and E is empty
        %warning('cmstep : Intersection of P and E are empty')
        x1  = [];
        L1  = [];
        bad = 2;
     else
        alpha = (A0'*x0 - B0)/A0_N2;
        x1    = x0 - ( (1+dim*alpha)/((dim+1)*A0_N2) )*A0hat;
        L1    = (dim^2*(1-alpha^2)/(dim^2-1))*(L0 - (2*(1+dim*alpha)/((dim+1)*(1+alpha)*A0_N))*AdM);
        bad   = 0;
     end

else
     A0hat   = L0*A0;
     A0_N    = A0'*A0hat;
     AdM     = A0hat*A0hat';
     x1      = x0 - (1/((dim+1)*sqrt(A0_N)))*A0hat;
     L1      = (dim^2/(dim^2-1))*(L0 - (2/((dim+1)*A0_N))*AdM);
     bad     = 0;
end

