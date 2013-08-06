function iqc_gui_test(filename)

if nargin==0,
    clc
    clear
    
    o_dir = pwd;
    
    test1_dir=which('CDelay_test.mdl');
    cd(test1_dir(1:end-15))
    
    all_file_name=dir;
    for i1=1:length(all_file_name)
        all_file{i1,1}=all_file_name(i1).name(1:end-4);
        if ~isempty(all_file{i1,1})
            test_file(all_file{i1,1})
        end
        if exist([all_file{i1,1},'.mdl.r14'])
            delete([all_file{i1,1},'.mdl.r14'])
        end
    end
    
    test1_dir=which('Sim_Deadzone_iqctest.mdl');
    cd(test1_dir(1:end-24))
    
    all_file_name=dir;
    for i1=1:length(all_file_name)
        all_file{i1,1}=all_file_name(i1).name(1:end-4);
        if ~isempty(all_file{i1,1})
            test_file(all_file{i1,1})
        end
        if exist([all_file{i1,1},'.mdl.r14'])
            delete([all_file{i1,1},'.mdl.r14'])
        end
    end
    
    cd(o_dir);
    clc
    clear
    disp('iqc_gui_test test done')
    
elseif nargin==1,
    test_file(filename)
else
    error('input error')
end

function test_file(filename)
clc
disp(filename)
filename_full=[filename '.mdl'];
open_system(filename_full)
sim(filename_full)
save_system(filename_full)
% disp('pause............................')
% pause
close_system(filename_full)
clear