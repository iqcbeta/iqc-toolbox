function y=m2s(x)
% function y=m2s(x)
%
% matrix to string conversion

[n,m]=size(x);
y='[';
for i=1:n,
    for j=1:m,
        y=[y  sprintf('%1.20g',x(i,j)) ' '];
    end
    y(length(y))=';';
end
y(length(y))=']';
