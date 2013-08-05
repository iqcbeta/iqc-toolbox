function [xc, opt, Del2F, V, count, eta] = center(A,b,xfeas)

% find the analytic center of a polyhedral
% defined by Ax <= b.

if nargin < 3 
   xfeas_avail = 0;
else
   xfeas_avail = 1;
end

% ---- find a feasible point ----
if ~xfeas_avail
   nA = size(A,1);
   mA = size(A,2);

   beta = 0.5;
   x = zeros(mA,1);
   pos = find(b<=1);

   bmod = b;
   bmod(pos) = 1;

   if ~isempty(pos) 
      done = 0;
   else 
      done = 1;
   end
      
   while ~done
     xa = ana(A,bmod,x);
     brep = max(b(pos), beta*bmod(pos) + (1-beta)*(A(pos,:)*xa));
     bmod(pos) = brep;
     x = xa;
     if all(bmod == b)
        done = 1;
     end
   end

   [xc,count,eta] = ana(A,b,x);
else
   [xc,count,eta] = ana(A,b,xfeas);
end

mA = size(A,1);
nA = size(A,2);
E = b-A*xc;
Del2F = inv(A'*diag(E.^(-2))*A);
V = (size(A,1)^size(A,2))/sqrt(det(A'*diag(E.^(-2))*A));
opt = -sum(log(E));



function [xa,count,eta] = ana(A,b,x) 

s = b-A*x;
if any(find(s<=0))
   s 
   error('x infesaible') 
end

done = 0;
count = 1;
   
while ~done
   sinv = 1./s;
   AS = zeros(size(A,2),length(sinv));
   for flcnt = 1:length(sinv)
       AS(:,flcnt) = A(flcnt,:)'*sinv(flcnt);
   end


   ND = -inv(AS*AS')*AS;
   dx = sum(ND,2);
   ds = -sum(A*ND,2);
   eta = sinv'*ds;

   if abs(eta)<1e-6    
      done = 1;
   else
      if eta<1
         s = s + ds;
         x = x + dx;
      else
         s = s + ds./(1+eta); 
         x = x + dx./(1+eta);
      end
   end
   count = count + 1;
end
xa = x;
