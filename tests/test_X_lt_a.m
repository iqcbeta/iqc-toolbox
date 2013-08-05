abst_init_lmi
t=rectangular;
[t 0;0 t]>1;
lmi_mincx_tbx(t)
disp(' ') 
disp(' The answer should be (very) close to 2; and there should be a warning message')
disp(' ')


abst_init_fdlmi
t=rectangular;
[t 0;0 t]>1;
fdlmi_mincx_tbx(t)
disp(' ') 
disp(' The answer should be (very) close to 2; and there should be a warning message')
disp(' ')
