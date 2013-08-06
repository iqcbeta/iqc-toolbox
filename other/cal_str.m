function out=cal_str(varargin)

out=0;
for i1=1:nargin
    out=out+str2double(varargin{i1});
end
out=num2str(out);