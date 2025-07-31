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
wn=3;
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

zpd=-imag(Sp(1))/(tan(D))+real(Sp(1));
PD=tf([1 -zpd],1);
figure(1); hold on;
rlocus(PD*G)
sgrid(zeta,wn)
%[Kpd poles]=rlocfind(PD*G)
Kpd=174.2554;
Gpd=feedback(Kpd*PD*G,1);
margin(Kpd*PD*G)
figure(2); hold on;
step(Gpd)
stepinfo(Gpd)
[y,t]=step(Gpd);
s=tf('s');
figure(3); hold on;
step(Gpd/s,10)
step(1/s,10)
[yr,tr]=step(Gpd/s,10);
nyquist(Kpd*PD*G)
