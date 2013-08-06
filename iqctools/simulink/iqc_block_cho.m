function iqc_block_cho(blk,block_name)

switch block_name
    %% L2signal
    case 'L2signal_type'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'sig_type')
            case 'Step'
                en{5} = 'on';
                en{6} = 'off';
                en{7} = 'off';
            case 'Square wave'
                en{5} = 'on';
                en{6} = 'on';
                en{7} = 'off';
            case 'Sine wave'
                en{5} = 'on';
                en{6} = 'on';
                en{7} = 'off';
            case 'Custom'
                en{5} = 'off';
                en{6} = 'off';
                en{7} = 'on';
            case 'Random'
                en{5} = 'on';
                en{6} = 'off';
                en{7} = 'off';
            case 'External'
                en{5} = 'off';
                en{6} = 'off';
                en{7} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
        %% performance
    case 'performance_imp'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'imp_met')
            case 'Only Simulation'
                en{2} = 'off';
                en{3} = 'off';
                en{4} = 'off';
            otherwise
                en{2} = 'on';
                en{3} = 'on';
                if strcmp(get_param(blk,'isiqcbode'),'on')
                    en{4} = 'on';
                end
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'performance_bodefre'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'isiqcbode')
            case 'on'
                en{4} = 'on';
            case 'off'
                en{4} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
        %% sim_saturation
    case 'sim_saturation_iqc1'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc1')
            case 'on'
                %                 en{5} = 'on';
                en{5} = 'on';
                en{6} = 'on';
                en{7} = 'on';
            case 'off'
                %                 en{5} = 'off';
                en{5} = 'off';
                en{6} = 'off';
                en{7} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_saturation_iqc2'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc2')
            case 'on'
                en{9} = 'on';
            case 'off'
                en{9} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_saturation_iqc3'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc3')
            case 'on'
                en{11} = 'on';
            case 'off'
                en{11} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_saturation_iqc4'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc4')
            case 'on'
                en{13} = 'on';
                en{14} = 'on';
            case 'off'
                en{13} = 'off';
                en{14} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_saturation_iqc5'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc5')
            case 'on'
                en{16} = 'on';
            case 'off'
                en{16} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
        %% sim_deadzone
    case 'sim_deadzone_iqc1'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc1')
            case 'on'
                %                 en{5} = 'on';
                en{5} = 'on';
                en{6} = 'on';
                en{7} = 'on';
            case 'off'
                %                 en{5} = 'off';
                en{5} = 'off';
                en{6} = 'off';
                en{7} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_deadzone_iqc2'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc2')
            case 'on'
                en{9} = 'on';
            case 'off'
                en{9} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_deadzone_iqc3'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc3')
            case 'on'
                en{11} = 'on';
            case 'off'
                en{11} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_deadzone_iqc4'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc4')
            case 'on'
                en{13} = 'on';
                en{14} = 'on';
            case 'off'
                en{13} = 'off';
                en{14} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_deadzone_iqc5'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc5')
            case 'on'
                en{16} = 'on';
            case 'off'
                en{16} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
        %% sim_intencdeadzone
    case 'sim_intencdeadzone_iqc1'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc1')
            case 'on'
                en{5} = 'on';
                en{6} = 'on';
                en{7} = 'on';
            case 'off'
                en{5} = 'off';
                en{6} = 'off';
                en{7} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_intencdeadzone_iqc2'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc2')
            case 'on'
                en{9} = 'on';
                en{10} = 'on';
                en{11} = 'on';
            case 'off'
                en{9} = 'off';
                en{10} = 'off';
                en{11} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_intencdeadzone_iqc3'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc3')
            case 'on'
                en{13} = 'on';
            case 'off'
                en{13} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
        %% sim_uncertaindelay
    case 'sim_uncertaindelay_sou'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'sou_del')
            case 'Random'
                en{3} = 'on';
                en{4} = 'off';
            case 'Custom'
                en{3} = 'off';
                en{4} = 'on';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_uncertaindelay_iqc1'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc1')
            case 'on'
                en{6} = 'on';
                en{7} = 'on';
            case 'off'
                en{6} = 'off';
                en{7} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_uncertaindelay_iqc2'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc2')
            case 'on'
                en{9} = 'on';
                en{10} = 'on';
            case 'off'
                en{9} = 'off';
                en{10} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_uncertaindelay_iqc3'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc3')
            case 'on'
                en{12} = 'on';
            case 'off'
                en{12} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
        %% sim_uncertaincdelay
    case 'sim_uncertaincdelay_sou'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'sou_del')
            case 'Random'
                en{3} = 'on';
                en{4} = 'off';
            case 'Custom'
                en{3} = 'off';
                en{4} = 'on';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_uncertaincdelay_iqc1'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc1')
            case 'on'
                en{6} = 'on';
                en{7} = 'on';
            case 'off'
                en{6} = 'off';
                en{7} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_uncertaincdelay_iqc2'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc2')
            case 'on'
                en{9} = 'on';
            case 'off'
                en{9} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
        %% sim_normboundsys
    case 'sim_normboundsys_sou'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'sou_sys')
            case 'Random'
                en{2} = 'on';
                en{3} = 'on';
                en{4} = 'on';
                en{5} = 'on';
                en{7} = 'off';
            case 'Custom'
                en{2} = 'off';
                en{3} = 'off';
                en{4} = 'off';
                en{5} = 'on';
                en{7} = 'on';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_normboundsys_iqc1'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc1')
            case 'on'
                en{10} = 'on';
                en{11} = 'on';
            case 'off'
                en{10} = 'off';
                en{11} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_normboundsys_iqc2'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc2')
            case 'on'
                en{13} = 'on';
            case 'off'
                en{13} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
        %% sim_constantgain
    case 'sim_constantgain_sou'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'sou_con')
            case 'Random'
                en{2} = 'on';
                en{3} = 'off';
            case 'Custom'
                en{2} = 'off';
                en{3} = 'on';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_constantgain_iqc1'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc1')
            case 'on'
                en{5} = 'on';
                en{6} = 'on';
            case 'off'
                en{5} = 'off';
                en{6} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_constantgain_iqc2'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc2')
            case 'on'
                en{8} = 'on';
            case 'off'
                en{8} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_constantgain_iqc3'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc3')
            case 'on'
                en{10} = 'on';
            case 'off'
                en{10} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_constantgain_iqc4'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc4')
            case 'on'
                en{12} = 'on';
            case 'off'
                en{12} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
        %% sim_ratelim
    case 'sim_ratelim_iqc1'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc1')
            case 'on'
                en{6} = 'on';
                en{7} = 'on';
            case 'off'
                en{6} = 'off';
                en{7} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_ratelim_iqc2'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc2')
            case 'on'
                en{9} = 'on';
            case 'off'
                en{9} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
        
        %% sim_window
    case 'sim_window_sou'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'sou_del')
            case 'Random'
                en{3} = 'on';
                en{4} = 'off';
            case 'Custom'
                en{3} = 'off';
                en{4} = 'on';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_window_iqc1'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc1')
            case 'on'
                en{7} = 'on';
                en{8} = 'on';
            case 'off'
                en{7} = 'off';
                en{8} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
    case 'sim_window_iqc2'
        en = get_param(blk,'MaskEnables');
        switch get_param(blk,'iqc2')
            case 'on'
                en{10} = 'on';
            case 'off'
                en{10} = 'off';
        end
        set_param(blk,'MaskEnables',en)
        return
end