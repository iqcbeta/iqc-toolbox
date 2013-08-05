function [w,xa,xb,xc,xd,dd]=iqc_d_slope(v,a,alpha,beta,ign)
% function [w,xa,xb,xc,xd,dd]=iqc_d_slope(v,a,alpha,beta,ign)
%
% defines iqc's for the relation w(t)=PHI(v(t)), 
% where PHI is a DIAGONAL operator:
%
%  PHI= [phi 0......0;
%        0  phi 0...0;
%        :          :
%        0......0 phi];    and phi is odd and satisfies:
%    
%           phi(x) - phi(y)      
% alpha  < ---------------- < alpha+beta     for beta>0 
%                x - y
%
% The IQC's have the form:
%
% \int q' D p  dt > \int q' H p dt
%
% where  p = (1+alpha/beta)v - 1/beta w   (auxiliar input)
%        q =  -alpha v +  w               (auxiliar output)
%
% D is diagonal and H=\H^T is a convolution operator 
% with:     sum_j ||h_ij||_1<=D_ii,   h_ij(t)>0 for all t
%
% v : input signal
% a : SYMMETRIC (n x n) cell of vectors a_ij  containing the (-)poles of h_ij
%
% IF a IS EITHER EMPTY OR NOT SPECIFIED, THEN H_ij=1 (static multipliers)
% IF a IS AN n x n  CELL, THEN 
%  (1) if a{i,j}(k)==Inf         H_ij(k) = 1                           
%  (2) if a{i,j}(k)==a           H_ij(k) = a/(s + a)   
%  (3) if a{i,j}(k)==b+i*c   i)  H_ij(k) = b/(s+b)+b*c/((s + b)^2+c^2)
%                            ii) H_ij(k+1) = b/(s+b)+b(s+b)/((s + b)^2+c^2) 
%
% ign=1 : ignores multipliers (3)-i
% ign=2 : ignores multipliers (3)-ii
% ign=0 : uses both multipliers (3)-i and (3)-ii
%
% default:   H_ij=1, alpha=0,  beta=1 ign=0 
%
% see also: IQC_D_SLOPE, IQC_SLOPE_ODD, IQC_SLOPE, IQC_MONOTONIC

%
% Written by fdamato@lids.mit.edu       June, 1998
% last modified Sept 19, 2000. July 24,1998 by fdamato@ecn.purdue,edu
%

n=size(v,1);

if nargin<5;ign=0;end
if ign<0|ign>2; error('allowed values for ign : 0,1,2');end
if nargin<4|isempty(beta);beta=1;end
if nargin<3|isempty(alpha);alpha=0;end
if nargin<2|isempty(a);
   a=cell(n,n);
   for n1=1:n;for n2=1:n;
     a{n1,n2}=Inf;
   end;end
else 
   if size(a,1)~=n|size(a,2)~=n;
      error('inadequate dimensions of argument "a"');
   end
   for n1=1:n;for n2=1:n;
      if (isempty(a{n1,n2})+isempty(a{n2,n1}))==0;
         if sort(a{n1,n2})~=sort(a{n2,n1});
            error('cell "a" must be symmetric');
         end 
      elseif (isempty(a{n1,n2})+isempty(a{n2,n1}))==1;
         error('cell "a" must be symmetric');
      end
   end;end  	
end
if nargin<1;error('input must be specified');end

s = tf([1 0],1);
w = signal(n);
r=1+alpha/beta;
p=r*v-1/beta*w;
q=-alpha*v+w;

Ninf=0;Nreal=0;Ncpx=0;
for n1=1:n;	% determine # of variables to define
   for n2=n1:n;  % changed here <----
      aux=length(a{n1,n2});
      for ndx=1:aux;
         if a{n1,n2}(ndx)==Inf;
            Ninf=Ninf+1;
         elseif isreal(a{n1,n2}(ndx)); 
            Nreal=Nreal+1;
         elseif ign==0
            Ncpx=Ncpx+2;
         else % ign~=0
            Ncpx=Ncpx+1;		
         end
      end % ndx
   end % n2
end % n1

% variables definitions and positiveness
dd=rectangular(1,n);		
xa=rectangular(1,Ninf+Nreal+Ncpx);
for ndx=1:size(xa,2);		
   xa{ndx}>0;
end

if (Nreal+Ncpx)>0;
   xc=rectangular(1,Nreal+Ncpx);
   for ndx=1:size(xc,2);
      xc{ndx}>0;
   end
end
if Ncpx>0; 
   xb=rectangular(1,Ncpx);
   xd=rectangular(1,Ncpx);
   for ndx=1:size(xb,2);
      xb{ndx}>0;
      xd{ndx}>0;
   end
end

countA=0;
countC=0;
countBD=0;

for n1=1:n
   n2=n1;  
   i=num2str(n1);
   q(n1)'*dd(n1)*p(n1)>0;
   eval(['sumROW_' i '=dd(n1);']);

      m=length(a{n1,n2}); % ---------diagonal entries 
      for ndx=1:m;
         aa=a{n1,n2}(ndx);
         countA=countA+1;
         if aa==Inf;		% static multiplier
            q(n1)'*-xa{countA}*p(n2)>0;
            eval(['sumROW_' i '=sumROW_' i '-xa{countA};']);
         elseif isreal(aa) % real pole
            if aa<0;
               aa=-aa;
               warning('element of "a" changed to be positive');
            end
            countC=countC+1;
            Hij=aa/(s+aa);          
            HijPj=Hij*p(n2);
            HijQj=Hij*q(n2);
      
            q(n1)'*-xa{countA}*HijPj > 0;
            p(n1)'*-xc{countC}*HijQj > 0;
            eval(['sumROW_' i '=sumROW_' i '-xa{countA};']);
            eval(['sumROW_' i '=sumROW_' i '-xc{countC};']);
         else 	% complex poles
            countC=countC+1;
            countBD=countBD+1;
            b=real(aa);b=abs(b);   %<<<== NOTE: b,c forced to be positive
            c=imag(aa);c=abs(c);
            s1_ij=b/(s+b);
            if ign~=1;
               s2_ij=b*c/((s+b)*(s+b)+c*c);
               Hij=s1_ij+s2_ij;
               Mij=s1_ij-s2_ij;
               NormHij=1+b*c/(b*b+c*c);
               NormMij=1-b*c/(b*b+c*c);
               HijPj=Hij*p(n2);
               HijQj=Hij*q(n2);
               MijPj=Mij*p(n2);
               MijQj=Mij*q(n2);

               q(n1)'*-xa{countA}*HijPj > 0;
               p(n1)'*-xc{countC}*HijQj > 0;
               eval(['sumROW_' i '=sumROW_' i '-xa{countA}*NormHij;']);
               eval(['sumROW_' i '=sumROW_' i '-xc{countC}*NormHij;']);
               q(n1)'*-xb{countBD}*MijPj > 0;
               p(n1)'*-xd{countBD}*MijQj > 0;
               eval(['sumROW_' i '=sumROW_' i '-xb{countBD}*NormMij;']);
               eval(['sumROW_' i '=sumROW_' i '-xd{countBD}*NormMij;']);
            end
            if ign==0;
               countA=countA+1;
               countC=countC+1;
               countBD=countBD+1;
            end
            if ign~=2;
               s3_ij=b*(s+b)/((s+b)*(s+b)+c*c);
               Rij=s1_ij+s3_ij;
               Tij=s1_ij-s3_ij;
               NormRij=1+b*b/(b*b+c*c);
               NormTij=1-b*b/(b*b+c*c);
               RijPj=Rij*p(n2);
               RijQj=Rij*q(n2);
               TijPj=Tij*p(n2);
               TijQj=Tij*q(n2);

               q(n1)'*-xa{countA}*RijPj > 0;
               p(n1)'*-xc{countC}*RijQj > 0;
               eval(['sumROW_' i '=sumROW_' i '-xa{countA}*NormRij;']);
               eval(['sumROW_' i '=sumROW_' i '-xc{countC}*NormRij;']);
               q(n1)'*-xb{countBD}*TijPj > 0;
               p(n1)'*-xd{countBD}*TijQj > 0;
               eval(['sumROW_' i '=sumROW_' i '-xb{countBD}*NormTij;']);
               eval(['sumROW_' i '=sumROW_' i '-xd{countBD}*NormTij;']);
            end	% ign     
         end % m 
      end % ndx

      
      for n2=n1+1:n
      m=length(a{n1,n2}); % --------- off diagonal entries 
      for ndx=1:m;
         aa=a{n1,n2}(ndx);
         countA=countA+1;
         if aa==Inf;		% static multiplier
            q(n1)'*-xa{countA}*p(n2)>0;
            q(n2)'*-xa{countA}*p(n1)>0;   % transpose counterpart
            eval(['sumROW_' i '=sumROW_' i '-xa{countA};']);
         elseif isreal(aa) % real pole
            if aa<0;
               aa=-aa;
               warning('element of "a" changed to be positive');
            end
            countC=countC+1;
            Hij=aa/(s+aa);          
            HijPj=Hij*p(n2);
            HijQj=Hij*q(n2);
      
            q(n1)'*-xa{countA}*HijPj > 0;
            p(n1)'*-xc{countC}*HijQj > 0;
	    
	    HjiPi=Hij*p(n1); % transpose counterparts
	    HjiQi=Hij*q(n1);
	    q(n2)'*-xa{countA}*HjiPi > 0;
            p(n2)'*-xc{countC}*HjiQi > 0;
	    	    
            eval(['sumROW_' i '=sumROW_' i '-xa{countA};']); % diag dominance
            eval(['sumROW_' i '=sumROW_' i '-xc{countC};']);
         else 	% complex poles
            countC=countC+1;
            countBD=countBD+1;
            b=real(aa);b=abs(b);   %<<<== NOTE: b,c forced to be positive
            c=imag(aa);c=abs(c);
            s1_ij=b/(s+b);
            if ign~=1;
               s2_ij=b*c/((s+b)*(s+b)+c*c);

               Hij=s1_ij+s2_ij;
               NormHij=1+b*c/(b*b+c*c);

               HijPj=Hij*p(n2);
               HijQj=Hij*q(n2);
               HjiPi=Hij*p(n1);
               HjiQi=Hij*q(n1);
	       
               q(n1)'*-xa{countA}*HijPj > 0;
               p(n1)'*-xc{countC}*HijQj > 0;
               q(n2)'*-xa{countA}*HjiPi > 0;
               p(n2)'*-xc{countC}*HjiQi > 0;
	       
	       
               eval(['sumROW_' i '=sumROW_' i '-xa{countA}*NormHij;']);
               eval(['sumROW_' i '=sumROW_' i '-xc{countC}*NormHij;']);

               Mij=s1_ij-s2_ij;
               NormMij=1-b*c/(b*b+c*c);
	       
               MijPj=Mij*p(n2);
               MijQj=Mij*q(n2);
               MjiPi=Mij*p(n1);
               MjiQi=Mij*q(n1);
	       
               q(n1)'*-xb{countBD}*MijPj > 0;
               p(n1)'*-xd{countBD}*MijQj > 0;
               q(n2)'*-xb{countBD}*MjiPi > 0;
               p(n2)'*-xd{countBD}*MjiQi > 0;
	       

               eval(['sumROW_' i '=sumROW_' i '-xb{countBD}*NormMij;']);
               eval(['sumROW_' i '=sumROW_' i '-xd{countBD}*NormMij;']);
            end
            if ign==0;
               countA=countA+1;
               countC=countC+1;
               countBD=countBD+1;
            end
            if ign~=2;
               s3_ij=b*(s+b)/((s+b)*(s+b)+c*c);
               Rij=s1_ij+s3_ij;
               NormRij=1+b*b/(b*b+c*c);	       

               RijPj=Rij*p(n2);
               RijQj=Rij*q(n2);	       
               RjiPi=Rij*p(n1);
               RjiQi=Rij*q(n1);	       

               q(n1)'*-xa{countA}*RijPj > 0;
               p(n1)'*-xc{countC}*RijQj > 0;
               q(n2)'*-xa{countA}*RjiPi > 0;
               p(n2)'*-xc{countC}*RjiQi > 0;
               eval(['sumROW_' i '=sumROW_' i '-xa{countA}*NormRij;']);
               eval(['sumROW_' i '=sumROW_' i '-xc{countC}*NormRij;']);
		     
               Tij=s1_ij-s3_ij;
               NormTij=1-b*b/(b*b+c*c);

               TijPj=Tij*p(n2);
               TijQj=Tij*q(n2);
               TjiPi=Tij*p(n1);
               TjiQi=Tij*q(n1);
	       
               q(n1)'*-xb{countBD}*TijPj > 0;
               p(n1)'*-xd{countBD}*TijQj > 0;
               q(n2)'*-xb{countBD}*TjiPi > 0;
               p(n2)'*-xd{countBD}*TjiQi > 0;	       
               eval(['sumROW_' i '=sumROW_' i '-xb{countBD}*NormTij;']);
               eval(['sumROW_' i '=sumROW_' i '-xd{countBD}*NormTij;']);
            end	% ign     
         end % m 
      end % ndx	
      end % n2


   for nfix1=1:n1-1  % contributions to sum_i from the lower triangular
     [countA,countC,countBD]=findex(nfix1,n1,a,ign);   
     aa = a{nfix1,n1};
     for ndx=1:length(aa);
        if aa(ndx)==Inf;
          countA=countA+1;	  
          eval(['sumROW_' i '=sumROW_' i '-xa{countA};']);
        elseif isreal(aa(ndx))
          countA=countA+1;	  
          countC=countC+1;
          eval(['sumROW_' i '=sumROW_' i '-xa{countA};']);
          eval(['sumROW_' i '=sumROW_' i '-xc{countC};']);
        elseif ~isreal(aa(ndx))
          countA=countA+1;	  
	  countC=countC+1;
          countBD=countBD+1;
	  b=real(aa);b=abs(b);   %<<<== NOTE: b,c forced to be positive
          c=imag(aa);c=abs(c);
	  
	  NormHij=1+b*c/(b*b+c*c);
          NormMij=1-b*c/(b*b+c*c);	  
          NormRij=1+b*b/(b*b+c*c);	       	  
          NormTij=1-b*b/(b*b+c*c);
	  
          if ign==0
             eval(['sumROW_' i '=sumROW_' i '-xa{countA}*NormHij;']);
             eval(['sumROW_' i '=sumROW_' i '-xc{countC}*NormHij;']);
             eval(['sumROW_' i '=sumROW_' i '-xb{countBD}*NormMij;']);
             eval(['sumROW_' i '=sumROW_' i '-xd{countBD}*NormMij;']);

	     countA=countA+1;
	     countC=countC+1;
	     countBD=countBD+1;
	     
             eval(['sumROW_' i '=sumROW_' i '-xa{countA}*NormRij;']);
             eval(['sumROW_' i '=sumROW_' i '-xc{countC}*NormRij;']);
             eval(['sumROW_' i '=sumROW_' i '-xb{countBD}*NormTij;']);
             eval(['sumROW_' i '=sumROW_' i '-xd{countBD}*NormTij;']);
 
          elseif ign==1
             eval(['sumROW_' i '=sumROW_' i '-xa{countA}*NormRij;']);
             eval(['sumROW_' i '=sumROW_' i '-xc{countC}*NormRij;']);
             eval(['sumROW_' i '=sumROW_' i '-xb{countBD}*NormTij;']);
             eval(['sumROW_' i '=sumROW_' i '-xd{countBD}*NormTij;']);
		   
          elseif ign==2
             eval(['sumROW_' i '=sumROW_' i '-xa{countA}*NormHij;']);
             eval(['sumROW_' i '=sumROW_' i '-xc{countC}*NormHij;']);
             eval(['sumROW_' i '=sumROW_' i '-xb{countBD}*NormMij;']);
             eval(['sumROW_' i '=sumROW_' i '-xd{countBD}*NormMij;']);	              
	  end % ign==0
        end % if aa==Inf
     end % ndx
   end % nfix
end % n1


for n1=1:n;
   i=num2str(n1);
   eval(['sumROW_' i '>0;']);
end

if (Nreal+Ncpx)==0; % assigns a(ny) value to xc to avoid warnings 
   xc=0;
end
if (Ncpx)==0; % assigns a(ny) value to xb, xd to avoid warnings 
   xb=0;xd=0;
end

return
%----- auxiliay subroutine --------

function  [countA,countC,countBD]=findex(ndx1,ndx2,a,ign);

if ndx1>ndx2;
  countAB=-1;
  countCD=-1;
  return
end
n=length(a);
countA=0;
countC=0;
countBD=0;

for n2=1:n
for n1=1:ndx1-1
  aa=a{n1,n2}
  m=length(aa);
  for ndx=1:m
    if aa(ndx)==Inf;
      countA=countA+1;
    elseif isreal(aa(ndx))
      countA=countA+1;
      countC=countC+1;
    elseif ~isreal(aa(ndx))
      if ign==0
        countA=countA+2;
        countC=countC+2;
	countBD=countBD+2;
      else 
        countA=countA+1;
        countC=countC+1;
	countBD=countBD+1;
      end
    end       
  end
end
end

for n2=ndx1:ndx2-1
  aa=a{ndx1,n2};
  m=length(aa);
  for ndx=1:m
    if aa(ndx)==Inf;
      countA=countA+1;
    elseif isreal(aa(ndx))
      countA=countA+1;
      countC=countC+1;
    elseif ~isreal(aa(ndx))
      if ign==0
        countA=countA+2;
        countC=countC+2;
	countBD=countBD+2;
      else 
        countA=countA+1;
        countC=countC+1;	
        countBD=countBD+1;
      end
    end  
  end
end

return


