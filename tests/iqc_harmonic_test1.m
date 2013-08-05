% testing iqc_harmonic on the system
% y''+y'+(a+b*cos(w0*t))y=f

clear all
N=100;            % # of steps in numerical integration
a=1;
b=2.5;
w0=3;
s=tf([1 0],1);

abst_init_iqc;
f=signal;
w=signal;   % w=cos(w0*t)*x
y=(1/(s*s+s+a))*(f-b*w);
w==iqc_harmonic(y,w0,[0.7 0.8 5]);
gain=iqc_gain_tbx(f,y)
iqc_bode
pause

disp(' now calculate the monodromy matrix in the associated ')
disp(' Hamiltonian system for a 10 percent smaller gain ...' )
disp(' ')

g=600/gain^2;
h=2*pi/(w0*N);
t=linspace(0,2*pi/w0,N);
M(:,:,1)=eye(4);
for i=1:N,
    k=a+b*cos(t(i)+h/2);
    M(:,:,i+1)=expm([0 1 0 0;-k -1 0 g;-1 0 0 k;0 0 -1 1]*h)*M(:,:,i);
end
[V,z]=eig(M(:,:,N+1));
z=diag(z);
abs_z=abs(z)
if all(abs(abs_z-1)>0.001),
   [y,ind]=sort(abs_z);
   U=V(:,ind(1:2));
   m=zeros(1,N);
   for i=1:N,
      xi=M(:,:,i)*U;
      m(i)=max(abs(real(eig(xi(3:4,1:2)/xi(1:2,1:2)))));
   end
   plot(t,m); grid;
end



