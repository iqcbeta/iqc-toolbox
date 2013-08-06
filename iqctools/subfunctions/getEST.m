str{sc}=['Y=[T11+T12*inv(T22)*T12'' -T12*inv(T22);',...
    '-inv(T22)*T12'' inv(T22)];'];
eval(str{sc});
sc=sc+1;

nx1=num2str(size(X,1));
str{sc}=['UVT=(eye(',nx1,')-X*Y);'];
eval(str{sc});
sc=sc+1;

str{sc}='[U,S,V]=svd(UVT);';
eval(str{sc});
sc=sc+1;

str{sc}='U=U*S;';
eval(str{sc});
sc=sc+1;

ny1=num2str(size(hatK,1));
ny2=num2str(size(hatM,1));
nx2=num2str(size(L,2));
str{sc}=['M_E=inv([U X*hatB;zeros(',ny2,',',ny1,') eye(',ny2,')])*',...
    '[hatK*inv(T1)-X*hatA*Y L;hatM*inv(T1) N]*',...
    'inv([V'' zeros(',ny1,',',nx2,');hatC*Y eye(',nx2,')]);'];
eval(str{sc});
sc=sc+1;

str{sc}=['AE=M_E(1:',ny1,',','1:',ny1,');'];
eval(str{sc});
sc=sc+1;

str{sc}=['BE=M_E(1:',ny1,',',cal_str('1',ny1),':end);'];
eval(str{sc});
sc=sc+1;

str{sc}=['CE=-M_E(',cal_str('1',ny1),':end,','1:',ny1,');'];
eval(str{sc});
sc=sc+1;

str{sc}=['DE=-M_E(1+',ny1,':end,',cal_str('1',ny1),':end);'];
eval(str{sc});
sc=sc+1;