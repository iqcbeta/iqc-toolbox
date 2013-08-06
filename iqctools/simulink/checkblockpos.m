function pos=checkblockpos(block_path,param_name)
% check block parameter position

global ext_iqc_setvalues

if ~ischar(block_path) || ~ischar(param_name)
    disp_str(9,'"block_path" and "param_name"','char')
end

pos=0;

if ~isfield(ext_iqc_setvalues,param_name)
    eval(['ext_iqc_setvalues.',param_name,'.idx=1;'])
    pos=1;
    return
else
    str=['ext_iqc_setvalues.',param_name];
    if ~isfield(eval(str),'idx')
        eval(['ext_iqc_setvalues.',param_name,'.idx=1;'])
        pos=1;
        return
    end
    eval(['pos=ext_iqc_setvalues.',param_name,'.idx;'])
    for i1=1:pos
        str=['ext_iqc_setvalues.',param_name,'.path{i1}'];
        if strcmp(eval(str),block_path)
            pos=i1;
            return
        end
    end
    pos=pos+1;
    eval(['ext_iqc_setvalues.',param_name,'.idx=pos;'])
    return
end
