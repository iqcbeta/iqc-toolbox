function out=IQC_par(n,iqc_name,parameter)
ns=num2str(n);
nx=length(parameter.name);
clear i j
for z=1:nx,
    switch (class(parameter.value{1,z}{1,1})),
        case 'char'
            str=[parameter.name{z},'=''',parameter.value{1,z}{1,1},''';'];
        case 'double'
            str=[parameter.name{z},'=',mat2str(parameter.value{1,z}{1,1}),';'];
        case 'cell'
            str=[parameter.name{z},'=parameter.value{1,z}{1,1};'];
    end
    eval(str);
end
switch iqc_name
    case 'ext_iqc_performance'
        out=['In=in',ns,';Out=out',ns];
        
        if uselmitool==1,
            if ~isempty(lmi_par)
                out={out; ['setlmioptions(''lmilab'',',...
                    mat2str(lmi_par),');']};
            else
                out={out; 'setlmioptions(''lmilab'');'};
            end
        elseif uselmitool==2,
            if ~isempty(lmi_par)
                str=[];
                for i=1:length(lmi_par)
                    if isnumeric(lmi_par{i})
                        str=[str,num2str(lmi_par{i}),','];
                    else
                        str=[str,'''',lmi_par{i},''','];
                    end
                end
                out={out; ['setlmioptions(''yalmip'',',str(1:end-1),');']};
            else
                out={out; 'setlmioptions(''yalmip'');'};
            end
        end
    case 'ext_iqc_ltvnorm'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_ltvnorm(in',ns,',',...
                mat2str(out_dim),',',mat2str(max_amp),')'];
        else
            out=['out',ns,'==iqc_ltvnorm(in',ns,',',...
                mat2str(out_dim),',',mat2str(max_amp),')'];
        end
    case 'ext_iqc_ratelimiter'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_ratelimiter(in',ns,',',...
                mat2str(pole),',', mat2str(k),')'];
        else
            out=['out',ns,'==iqc_ratelimiter(in',ns,',',...
                mat2str(pole),',', mat2str(k),')'];
        end
    case 'ext_iqc_tvscalar'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_tvscalar(in',ns,',',...
                mat2str(upp_bou),')'];
        else
            out=['out',ns,'==iqc_tvscalar(in',ns,',',...
                mat2str(upp_bou),')'];
        end
    case 'ext_iqc_harmonic'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_harmonic(in',ns,',',...
                mat2str(del_fre),',',mat2str(pole),')'];
        else
            out=['out',ns,'==iqc_harmonic(in',ns,',',...
                mat2str(del_fre),',',mat2str(pole),')'];
        end
    case 'ext_iqc_sector'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_sector(in',ns,',',...
                mat2str(low_sec),',',mat2str(upp_sec),')'];
        else
            out=['out',ns,'==iqc_sector(in',ns,',',...
                mat2str(low_sec),',',mat2str(upp_sec),')'];
        end
    case 'ext_iqc_popov'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_popov(in',ns,',',...
                mat2str(low_sec),',',mat2str(upp_sec),')'];
        else
            out=['out',ns,'==iqc_popov(in',ns,',',...
                mat2str(low_sec),',',mat2str(upp_sec),')'];
        end
    case 'ext_iqc_popov_vect'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_popov_vect(in',ns,',''',...
                sign_par,''')'];
        else
            out=['out',ns,'==iqc_popov_vect(in',ns,',''',...
                sign_par,''')'];
        end
    case 'ext_iqc_monotonic'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_monotonic(in',ns,',',...
                mat2str(pole),',',mat2str(ope),',',mat2str(int_dd),')'];
        else
            out=['out',ns,'==iqc_monotonic(in',ns,',',...
                mat2str(pole),',',mat2str(ope),',',mat2str(int_dd),')'];
        end
    case 'ext_iqc_slowtv'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_slowtv(in',ns,',',...
                mat2str(upp_bou_dd),',',mat2str(pole),',',mat2str(upp_bou_d),')'];
        else
            out=['out',ns,'==iqc_slowtv(in',ns,',',...
                mat2str(upp_bou_dd),',',mat2str(pole),',',mat2str(upp_bou_d),')'];
        end
    case 'ext_iqc_window'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_window(in',ns,',',...
                mat2str(max_del),',',mat2str(pole),')'];
        else
            out=['out',ns,'==iqc_window(in',ns,',',...
                mat2str(max_del),',',mat2str(pole),')'];
        end
    case 'ext_iqc_polytope'
        nx=size(arr_pol,2);
        str=[];
        for i=1:nx
            str=[str,mat2str(arr_pol{i}),','];
        end
        str=['{',str(1:end-1),'}'];
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_polytope(in',ns,',',...
                str,')'];
        else
            out=['out',ns,'==iqc_polytope(in',ns,',',...
                str,')'];
        end
    case 'ext_iqc_polytope_stvp'
        nx=size(arr_del,2);
        str1=[];
        for i=1:nx
            str1=[str1,mat2str(arr_del{i}),','];
        end
        str1=['{',str1(1:end-1),'}'];
        nx=size(arr_ddel,2);
        str2=[];
        for i=1:nx
            str2=[str2,mat2str(arr_ddel{i}),','];
        end
        str2=['{',str2(1:end-1),'}'];
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_polytope_stvp(in',ns,',',...
                str1,',',str2,',',mat2str(str_mat),')'];
        else
            out=['out',ns,'==iqc_polytope_stvp(in',ns,',',...
                str1,',',str2,',',mat2str(str_mat),')'];
        end
    case 'ext_iqc_diag'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_diag(in',ns,');'];
        else
            out=['out',ns,'==iqc_diag(in',ns,')'];
        end
    case 'ext_iqc_ltiunmod'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_ltiunmod(in',ns,',',...
                mat2str(pole),',',mat2str(out_dim),',',mat2str(upp_bou),')'];
        else
            out=['out',ns,'==iqc_ltiunmod(in',ns,',',...
                mat2str(pole),',',mat2str(out_dim),',',mat2str(upp_bou),')'];
        end
    case 'ext_iqc_slope_odd'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_slope_odd(in',ns,',',...
                mat2str(pole),',',mat2str(len_def_h),',',...
                mat2str(low_bou),',',mat2str(upp_bou),')'];
        else
            out=['out',ns,'==iqc_slope_odd(in',ns,',',...
                mat2str(pole),',',mat2str(len_def_h),',',...
                mat2str(low_bou),',',mat2str(upp_bou),')'];
        end
    case 'ext_iqc_slope'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_slope(in',ns,',',...
                mat2str(pole),',',mat2str(len_def_h),',',...
                mat2str(low_bou),',',mat2str(upp_bou),')'];
        else
            out=['out',ns,'==iqc_slope(in',ns,',',...
                mat2str(pole),',',mat2str(len_def_h),',',...
                mat2str(low_bou),',',mat2str(upp_bou),')'];
        end
    case 'ext_iqc_dzn_e_odd'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_dzn_e_odd(in',ns,',',...
                mat2str(pole),',',mat2str(def_num_h),',',mat2str(slo_bou),')'];
        else
            out=['out',ns,'==iqc_dzn_e_odd(in',ns,',',...
                mat2str(pole),',',mat2str(def_num_h),',',mat2str(slo_bou),')'];
        end
    case 'ext_iqc_dzn_e'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_dzn_e(in',ns,',',...
                mat2str(pole),',',mat2str(def_num_h),',',mat2str(slo_bou),')'];
        else
            out=['out',ns,'==iqc_dzn_e(in',ns,',',...
                mat2str(pole),',',mat2str(def_num_h),',',mat2str(slo_bou),')'];
        end
    case 'ext_iqc_cdelay'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_cdelay(in',ns,',',...
                mat2str(max_tim),',',mat2str(pole),')'];
        else
            out=['out',ns,'==iqc_cdelay(in',ns,',',...
                mat2str(max_tim),',',mat2str(pole),')'];
        end
    case 'ext_iqc_delay'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_delay(in',ns,',',...
                mat2str(max_tim),',',mat2str(pole),')'];
        else
            out=['out',ns,'==iqc_delay(in',ns,',',...
                mat2str(max_tim),',',mat2str(pole),')'];
        end
    case 'ext_iqc_delay1'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_delay1(in',ns,',',...
                mat2str(max_tim),',',mat2str(pole),')'];
        else
            out=['out',ns,'==iqc_delay1(in',ns,',',...
                mat2str(max_tim),',',mat2str(pole),')'];
        end
    case 'ext_iqc_white'
        if ~isempty(outvariable)
            out=['[out',ns,',',outvariable(1:length(outvariable)-1),...
                ']=iqc_white(',num2str(sig_size(2)),',',...
                mat2str(ban_whi),',',mat2str(pole),')'];
        else
            out=['out',ns,'=iqc_white(',num2str(sig_size(2)),',',...
                mat2str(ban_whi),',',mat2str(pole),')'];
        end
    case 'ext_iqc_ltigain'
        try
            if ~isempty(addterm_times)
                if ~isempty(outvariable)
                    out=['[out',ns,',',...
                        outvariable(1:length(outvariable)-1),...
                        ']=iqc_ltigain(in',ns,',',...
                        mat2str(pole),')|',num2str(addterm_times),'*'];
                else
                    out=['out',ns,'==iqc_ltigain(in',ns,',',...
                        mat2str(pole),')|',num2str(addterm_times),'*'];
                end
            end
        catch err
            if ~isempty(outvariable)
                out=['[out',ns,',',...
                    outvariable(1:length(outvariable)-1),...
                    ']=iqc_ltigain(in',ns,',',...
                    mat2str(pole),')'];
            else
                out=['out',ns,'==iqc_ltigain(in',ns,',',...
                    mat2str(pole),')'];
            end
        end
    otherwise
        out='';
end
