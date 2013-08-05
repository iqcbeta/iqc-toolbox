setlmis([]);
[nvar,nlmivar,svar]=lmivar(2,[1 1]);
[nvar,nlmivar,svar]=lmivar(2,[1 1]);
[nvar,nlmivar,svar]=lmivar(2,[1 1]);
[nvar,nlmivar,svar]=lmivar(2,[1 1]);
[nvar,nlmivar,svar]=lmivar(3,1);
[nvar,nlmivar,svar]=lmivar(3,2);
[nvar,nlmivar,svar]=lmivar(3,3);
[nvar,nlmivar,svar]=lmivar(3,4);
[nvar,nlmivar,svar]=lmivar(3,1);
[nvar,nlmivar,svar]=lmivar(3,2);
[nvar,nlmivar,svar]=lmivar(3,3);
[nvar,nlmivar,svar]=lmivar(3,4);
[nvar,nlmivar,svar]=lmivar(1,[1 1]);
[nvar,nlmivar,svar]=lmivar(1,[2 1]);
lmiterm([1 1 1 9],[-1],[1],'s');
lmiterm([2 1 1 10],[-1],[1],'s');
lmiterm([3 1 1 11],[-1],[1],'s');
lmiterm([4 1 1 12],[-1],[1],'s');
lmiterm([5 1 1 14],[-1 0;0 -1],[1 0;0 1],'s');

pp=lmivar(1,[8 1]);
lmiterm([6 1 1 15],[eye(8);zeros(2,8)], ...
   [-50 0 0 0 0 0 0 0 0 8; 3.125 -2 -2.75 0 0 0 0 0 0.5 0; ...
      0 4 0 0 0 0 0 0 0 0; 0 0 128 -100.01 -0.125 0 0 0 0 0; ...
      0 0 0 8 0 0 0 0 0 0; 0 0 0 0 0 -50 0 0 0 8; ...
      0 0 -620 468.359375 -4.296875 0 -1 0 0 -1; ...
      0 0 0 0 0 -6.25 0 -1 0 1],'s');
lmiterm([6 1 1 -5], ...
 [0;0;-620;468.359375;-4.296874999999999;0;-1;0;0;-1], ...
 [0 0 0 0 0 -6.25 0 0 0 1],'s');
lmiterm([6 1 1 -6], ...
 [0;0;-620;468.359375;-4.296874999999999;0;1;0;0;-1], ...
 [0 0 0 0 0 -6.25 0 0 0 1],'s');
lmiterm([6 1 1 -7],[0;0;0;0;0;-6.25;0;-1;0;1], ...
 [0 0 -620 468.359375 -4.296874999999999 0 0 0 0 -1],'s');
lmiterm([6 1 1 -8],[0;0;0;0;0;-6.25;0;1;0;1], ...
 [0 0 -620 468.359375 -4.296874999999999 0 0 0 0 -1],'s');
lmiterm([6 1 1 -13], ...
 [0;-2480;59950;-46874.99609374999;-58.544921875;312.5;0;0;0;-50], ...
 [0 0 0 0 0 -6.25 0 0 0 1],'s');
lmiterm([6 1 1 -14], ...
   [0 0;0 0;-1240 0;936.71875 0;-8.593749999999998 0; ...
      0 6.25;0 0;0 0;0 0;-1 0], ...
   [0 0 -1240 936.71875 -8.593749999999998 0 0 0 0 -1; ...
      0 0 0 0 0 6.25 0 0 0 0],'s');
lmiterm([6 1 1 -14], ...
   [0 0;0 0;-1240 0;936.71875 0;-8.593749999999998 0; ...
      0 0;0 0;0 0;0 0;0 0], ...
   [0 0 620 -468.359375 4.296874999999999 0 0 0 0 2; ...
      0 0 0 0 0 0 0 0 0 0],'s');
a=zeros(10); a(3,3)=0.25;
lmiterm([6 1 1 0],a);
[gai,ndec,svar]=lmivar(1,[1 1]);
a=zeros(10); a(9,9)=1;
lmiterm([-6 1 1 16],a,1);

% commenting out this "lmiterm" makes the LMI feasible !!
% note that LMI's #5 and 7 are defined by single lmiterms,
% and are identical !!!
% What makes me mad, is that changing the location of
% the line
% lmiterm([5 1 1 14],[-1 0;0 -1],[1 0;0 1],'s');
% (moving it to this place from where it is now)
% will make the LMI's infeasible AGAIN !!!
% i.e. feasibility depends on the ORDER in which the
% LMIs are declared
% Also, changing the feasibility radius parameter to -1 (infinite)
% makes the system infeasible!
lmiterm([7 1 1 14],[-1 0;0 -1],[1 0;0 1],'s');

lmi=getlmis;
c=[zeros(1,44) 1];
[gain,xopt]=mincx(lmi,c,[0.0001 500 0 50 0]);
