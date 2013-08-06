function y=decs2mats(x,s)
% function y=decs2mats(x,s)
% x is a 1xN matrix, s is the structure matrix in the LMI Toolbox sense:
%       its entries are integers -N<=s(i,j)<=N;
% y is the matrix of same size as s, y(i,j)=sign(s(i,j))*x(abs(s(i,j))

if size(s,1)==1
    x=reshape(x,1,length(x));
end

y=sign(s).*x(abs(s)+(s==0));
