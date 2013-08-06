function varargout=iqctool(varargin)

switch varargin{1}
    case {'ver','version'}
        varargout{1}='iqcb036';
        
    case 'clear'
        clear global ABST
        clear global ABSTSOLUTION
        clear global ext_iqc_setvalues
        
    case 'simclear'
        clear global ext_iqc_setvalues
        
end