% default values

G=tf([-1 1],[1 2 5]);

G=G/(norm(G,inf)+1);

% calculating the lower bound for the gain
[m,p,ww]=bode(G);
ww=linspace(0,max(ww)/3,1000);
g=freqresp(G,ww);
g=squeeze(g);
lbw=abs(g./(1+g)).*(real(g)>-1)+abs(g./imag(g)).*(real(g)<=-1);
lb=max(lbw);


abst_init_iqc

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sedumi')
% setlmioptions('yalmip','solver','sdpt3')

f=signal;
w=signal;
v=G*(f-w);
w==iqc_popov(v,[0 1]);

setlmioptions('lmilab')
lmilab_v=iqc_gain_tbx(f,w);
data_v{1,1}='lmilab';
data_v{2,1}=lmilab_v;

if exist('sdpt3')~=0
    setlmioptions('yalmip','solver','sdpt3')
    sdpt3_v=iqc_gain_tbx(f,w);
    data_v{1,2}='sdpt3';
    data_v{2,2}=sdpt3_v;
end

if exist('sedumi')~=0
    setlmioptions('yalmip','solver','sedumi')
    sedumi_v=iqc_gain_tbx(f,w);
    data_v{1,3}='sedumi';
    data_v{2,3}=sedumi_v;
end

fprintf('\n')
fprintf('\n')
fprintf('----------------------------------------\n')
for i1=1:size(data_v,2)
    disp([data_v{1,i1},': ',num2str(data_v{2,i1})])
end
