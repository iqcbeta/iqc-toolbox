function varargout=value_iqc(varargin)
% function varargout=value_iqc(varargin)
%
% To retrieve optimal values of "abst" elements (varargin)
% and place them in varargout
%
% Written by ameg@mit.edu,  last modified October 13, 1997
%
% change name to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu
% Last modified by cmj on 2013/4/18

global ABSTSOLUTION
if isempty(ABSTSOLUTION),
    disp_str(13)
end

% if nargin~=nargout,
%    error('number of inputs must match number of outputs')
% end

for k=1:nargin,
    varargout{k}=ABSTSOLUTION{double(varargin{k})}; %#ok<*AGROW>
end
