  if nargin < 6
     error('Not enough input arguments ...')
  end

  x0  = varargin{1};
  L0  = varargin{2};
  A0  = varargin{3};
  B0  = varargin{4};
  Ind = varargin{5};
  dim = length(x0);
  Volorg = pi^(dim/2)*sqrt(abs(det(L0)))/gamma(dim/2+1);

  for flcnt1 = 1: length(Ind)
      g   = A0(Ind(flcnt1),:);
      gb  = B0(Ind(flcnt1));
      [x, L, badflag] = e_cut(x0,L0,[g';gb]);
      E_x{flcnt1,1}   = x;
      E_L{flcnt1,1}   = L; 
      if badflag == 0
         detL = det(L);
         if detL > 0
            E_Vol(flcnt1,1) = pi^(dim/2)*sqrt(detL)/gamma(dim/2+1);
         else 
            E_Vol(flcnt1,1) = 0;
         end
      elseif badflag == 1
         E_Vol(flcnt1,1) = Volorg;
      elseif badflag == 2
         E_Vol(flcnt1,1) = 0;
      end
  end

  ind1 = find(E_Vol > 0);
  E_Volr = E_Vol(ind1);
 
  if ~isempty(ind1)
     [minV, ind2] = min(E_Volr);
     x0   = E_x{ind2,1}; 
     L0   = E_L{ind2,1};
     V0   = minV;
     done = 0;
  else
     x0   = [];
     L0   = [];
     V0   = 0; 
     done = 1; 
  end

  while ~done
    cor = find(A0*x0-B0 >= 0);
    if ~isempty(cor)
        num = length(cor);
        E_x{num,1}  = [];
        E_L{num,1}  = [];
        E_Vol       = zeros(num,1); 
        for flcnt1 = 1: length(cor)
            g   = A0(cor(flcnt1),:);
            gb  = B0(cor(flcnt1));
            [x, L, badflag] = e_cut(x0,L0,[g';gb]);
            E_x{flcnt1,1}   = x;
            E_L{flcnt1,1}   = L; 
            if badflag == 0
               detL = det(L);
               if detL > 0
                  E_Vol(flcnt1,1) = pi^(dim/2)*sqrt(detL)/gamma(dim/2+1);
               else 
                  E_Vol(flcnt1,1) = 0;
               end
            elseif badflag == 1
               E_Vol(flcnt1,1) = Volorg;
            elseif badflag == 2
               E_Vol(flcnt1,1) = 0;
            end
        end

        ind1 = find(E_Vol > 0);
        E_Volr = E_Vol(ind1);
 
        if ~isempty(ind1)
           [minV, ind2] = min(E_Volr);
           x1   = E_x{ind2,1}; 
           L1   = E_L{ind2,1};
           V1   = minV;
        else
           x1   = [];
           L1   = [];
           V1   = 0; 
        end

	    if isempty(x1)
           done = 1;           
        else 
           if norm(x1 - x0) < 1e-3
              done = 1;
           else
              x0 = x1;
              L0 = L1;
              V0 = V1;
           end
        end    
    else    
        done = 1;
    end
  end

  x1 = x0;
  varargout(1) = {L0};
  varargout(2) = {V0};

end
