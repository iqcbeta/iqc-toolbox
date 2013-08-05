%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                      %%
%%  Hamiltonian system has eigenvalues  %%
%%  on the imaginary axis               %%
%%                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function omega=IS_Hsys_OK(A,B,Q,F,R)


H11=A-B*inv(R)*F';
H12=B*inv(R)*B';
H21=Q-F*inv(R)*F';
H22=-A'+F*inv(R)*B';
H=[H11,H12;H21,H22]; %'

DH      = eig(H);
pos     = find(abs(real(DH))<1e-6);
imag_DH = imag(DH);           % imaginery part of eigenvalues 
                              % which are on the imaginery axis

if ~isempty(pos)
    omega = imag_DH(pos);
    pos   = find(imag_DH(pos)>=0);
    omega = omega(pos);
else
    omega=[];
end



