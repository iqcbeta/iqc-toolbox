function setblockparam(pos,glparam_name,alparam_name,alparam) %#ok<*INUSL>

global ext_iqc_setvalues %#ok<*NUSED>

param_n=length(alparam);

for i1=1:param_n
    param_name=alparam_name{i1};
    str=['ext_iqc_setvalues.',glparam_name,'.',param_name,'{pos}=alparam{i1};'];
    eval(str);    
end
