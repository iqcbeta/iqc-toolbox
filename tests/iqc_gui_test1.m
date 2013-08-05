%
% This is a test program for GUI
% 
abst_init_iqc
lmitbx_options([0 0 0 0 1]);

disp('*** GUI test1 ...')
open_system('test1')
gain=iqc_gui('test1')
ga=1;
if isempty(gain)|(abs(gain-ga)>ga/1000),
   error('iqc_gui is corrupted')
end
close_system('test1')

disp('*** GUI test2 ...')
open_system('test2')
gain=iqc_gui('test2')
ga=1.8474;
if isempty(gain)|(abs(gain-ga)>ga/1000),
   error('iqc_gui is corrupted')
end
close_system('test2')

disp('*** GUI test3 ...')
open_system('test3')
gain=iqc_gui('test3')
ga=0.6903;
if isempty(gain)|(abs(gain-ga)>ga/1000),
   error('iqc_gui is corrupted')
end
close_system('test3')

disp('*** GUI test4 ...')
open_system('test4')
gain=iqc_gui('test4')
ga=0.1429;
if isempty(gain)|(abs(gain-ga)>ga/1000),
   error('iqc_gui is corrupted')
end
close_system('test4')

disp(' ')
disp(' ')
disp('***********iqc_gui.m TESTING OK!!!')
