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
num=[K*n];
numv=[K*n 0];
den=[L*Jo,L*bo+R*Jo,R*bo+K^2,0];
G=tf(num,den);
Po=4;
Ts=2;
zeta=0.75;
wn=4;
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

zL=-P(3);
PoleAngle=atan2(imag(Sp(1)+zL),real(Sp(1)+zL))-D;
pL=(imag(Sp(1))/tan(PoleAngle))-real(Sp(1));
Gc=tf([1 zL],[1 pL]);
T1=1/zL;
gamma=pL*T1;
rlocus(Gc*G)
sgrid(zeta,wn);
%[Kc,pole]=rlocfind(Gc*G);
Kc=1534.7 ;
Gc3=feedback(Kc*Gc*G,1);
step(Gc3)
stepinfo(Gc3)
[y,t]=step(Gc3);
s=tf('s');
figure(3); hold on;
step(Gc3/s)
step(1/s)
axis([0 7 0 7])
[yr,tr]=step(Gc3/s,7);
margin(Kc*Gc*G)
nyquist(Kc*Gc*G)