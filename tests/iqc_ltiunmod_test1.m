% testing "iqc_ltiunmod": iqc description for
% a LTI unmodeled dynamics with bounded norm
%
% The L2 gain from f to v in the interconnection 
%
% v = G0(f+w),  where v is 3x1, f and w are scalar
% signals, and w = sum_{over i} delta(i)*v(i)
% || delta || < k 
%
% is computed to check the stability of the uncertain 
% system. Also, mu-toolbox is used to verify the correctness 
% of the stability prediction. 
%
% Written by C.Kao cykao@mit.edu,  last modified July 14, 1999

disp(' ')
disp(' Sys : v=G0*w, w=sum_{over i} delta(i)*v(i) ')
disp(' G0  : A=[0.9501,0.4860,0.4565; 0.2311,0.8913,0.0185; 0.6068,0.7621,0.8214]')
disp('       B=[0.889;  1.231;  1.584]')
disp('       C=sys(3)')
disp('       D=[0;0;0]')
disp(' ')
disp(' First, use mu-toolbox to get maximun ||delta|| such that the system')
disp(' is still stable')

clear all
A=[0.9501, 0.4860, 0.4565;
   0.2311, 0.8913, 0.0185;
   0.6068, 0.7621, 0.8214];
B=[0.889;  1.231;  1.584];
C=eye(3);
D=[0;0;0];
G0=ss(A,B,C,D);

M=pck(A,B,eye(3),zeros(3,1));
omega=logspace(-2,2,100);
Mg=frsp(M,omega);
delta=[1 3];
mubnd=mu(Mg,delta,'lus');
k0=(1/max(mubnd(:,1)));
disp(' ')
disp(' mu-toolbox predicts that the ststem is stable when ')
disp([' the uncertainty has norm less than ',num2str(k0)])

abst_init_iqc
% lmitbx_options([0 0 0 0 1]);

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sdpt3')

k1=0.999*k0;
f=signal;
w=signal;
y=G0*(f+w);
w==iqc_ltiunmod(y,1,1,k1);
g1=iqc_gain_tbx(f,y);

abst_init_iqc
k2=0.9999*k0;
f=signal;
w=signal;
y=G0*(f+w);
w==iqc_ltiunmod(y,1,1,k2);
g2=iqc_gain_tbx(f,y);

disp(' ')
if ~isempty(g2) && ~isinf(g2)
   str1=' IQC theory also proves that system can substain an uncertainty with ';
   str2=' norm less than ';
   disp(str1)
   disp([str2,num2str(k2)])
elseif ~isempty(g1) && ~isinf(g1)
   str1=' IQC theory also proves that system can substain an uncertainty with ';
   str2=' norm less than ';
   disp(str1)
   disp([str2,num2str(k1)])
else
   str1=' IQC theory is NOT able to prove that the system is stable, ';
   str2=[' when the uncertainty has norm equal to ',num2str(k1)];
   disp(str1)
   disp(str2)
end


