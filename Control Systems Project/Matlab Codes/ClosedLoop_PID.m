close all; clear; clc;
Jm=0.04;
Jl=0.5;
bm=0.1;
bl=0;
K=0.25;
R=8;
L=0.05;
n=1/75;
Jo=Jm+(n^2)*Jl;
bo=bm+(n^2)*bl;
num=[K*n];;
numv=[K*n 0];
den=[L*Jo,L*bo+R*Jo,R*bo+K^2,0];
G=tf(num,den);
Po=4;
Ts=2;
[zeta,wn]=SecondOrderResponse(Po,Ts);
P=pole(G);
Z=zero(G);
%D=-pi+sum(P)-sum(Z)
if zeta<1
    Sp(1)=-zeta*wn + wn*sqrt(1-zeta^2)*j;
    Sp(2)=-zeta*wn - wn*sqrt(1-zeta^2)*j;
else
    Sp(1)=-zeta*wn + wn*sqrt(zeta^2-1)*j;
    Sp(2)=-zeta*wn - wn*sqrt(zeta^2-1)*j;
end
D=-pi;

for i=1:length(P)
    D=D+angle(Sp(1)-P(i));
end

for i=1:length(Z)
    D=D-angle(Sp(1)-Z(i));
end

zpid1=-5;
Dnew=D+angle(Sp(1))-angle(Sp(1)-zpid1)
zpi=-(imag(Sp(1))/tan(Dnew))+real(Sp(1))
PID=tf(conv([1 -zpid1],[1 -zpi]),[1 0])
figure(4);hold on;
rlocus(PID*G)
sgrid(zeta,wn)
%[Kpid poles]=rlocfind(PID*G)
Kpid=169.8392;
Gc2=feedback(Kpid*PID*G,1);
figure(1); hold on;
step(Gc2)
[y,t]=step(Gc2);
stepinfo(Gc2)
figure(2); hold on;
margin(Kpid*PID*G)
s=tf([1 0],[1]);
figure(3);hold on;
step(1/s)
step(Gc2/s)%steady state error is zero for ramp
axis([0 6 0 6])
[y,t]=step(Gc2/s);%ramp;
s2=tf([1 0 0],[1]);
figure(5);hold on;
step(2/s2)
step(2*Gc2/s2)%steady state error is zero for ramp
axis([0 12 0 144])
[y,t]=step(2*Gc2/(s2));%parabola has error, makes sense since type 2
nyquist(Kpid*PID*G)
