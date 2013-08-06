function run_lmi_exe
% function run_lmi_exe
%
% runs the LMI Control Toolbox code stored in lmi_exe.mat

load lmi_exe
for k=1:length(exe),
    disp([exe{k} ';'])
    eval([exe{k} ';'])
end
