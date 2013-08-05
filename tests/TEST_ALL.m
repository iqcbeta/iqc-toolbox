% TEST_ALL.m
%
% Purpose: perform all tests included in iqc_toolbox
%          There is a 'pause' somewhere in one of the tests but otherwise
%          will run automatically
%
% except for "manual_test_GS.m" -------> looks for "iqc_gain_tbx_YALMIP"
%            "iqc_d_slope_odd_test"
%            "iqc_d_slope_odd_test1"
%            "iqc_d_slope_odd_test2" ---> all three above gives following problem:
%                                        ??? Error using ==> iqc_extract
%                                        illegal differentiation log #115
%
% Andres Marcos, 27-October-2006, marcosa@ame.umn.edu

fdlmi_test
close all; clear all

iqc_beta_test
close all; clear all

iqc_cdelay_test
close all; clear all

iqc_delay1_test
close all; clear all

iqc_delay_test
close all; clear all

if exist('butter')~=0
    iqc_domharmonic_test
    close all; clear all
end

iqc_dzn_e_odd_test
close all; clear all

iqc_dzn_e_test
close all; clear all

iqc_gain_tbx_test1
close all; clear all

iqc_get_mlmi_test
close all; clear all

iqc_gui_test1
close all; clear all

iqc_harmonic_test1
close all; clear all

iqc_ltigain_test1
close all; clear all

iqc_ltigain_test2
close all; clear all

iqc_ltiunmod_test
close all; clear all

iqc_ltiunmod_test1
close all; clear all

iqc_ltvnorm_test1
close all; clear all

iqc_monotonic_test1
close all; clear all

% iqc_multi_harmonic_test
% close all; clear all

iqc_polytope_stvp_test
close all; clear all

iqc_polytope_stvp_test2
close all; clear all

iqc_polytope_test
close all; clear all

iqc_popov_test1
close all; clear all

iqc_popov_vect_test2
close all; clear all

iqc_sector_popov_vect_test
close all; clear all

iqc_slope_odd_test
close all; clear all

iqc_slope_test
close all; clear all

iqc_slowtv_test
close all; clear all

iqc_tvscalar_test
close all; clear all

iqc_value_test1
close all; clear all

iqc_white_test
close all; clear all

iqc_window_test
close all; clear all

lmi_feas_tbx_test1
close all; clear all

lmi_feas_tbx_test2
close all; clear all

lmi_mincx_tbx_test1
close all; clear all

lmi_mincx_tbx_test2
close all; clear all

lmitest_s
close all; clear all

manual_test1
close all; clear all

test_X_lt_a
close all; clear all



