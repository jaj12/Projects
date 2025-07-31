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
wn=2.79;
zeta=0.71;
numv=[K*n 0];
den=[L*Jo,L*bo+R*Jo,R*bo+K^2,0];
G=tf(num,den);
%%
K=628.8254;
Kc=1392.747287;
zero=1.747;
pole=3.869458;
Gc=tf([1 1.747],[1 3.869458]);
Gclosed=feedback(Kc*Gc*G,1);
figure(1); hold on;
margin(Kc*Gc*G)
figure(2);hold on;
step(Gclosed)
stepinfo(Gclosed)
[y,t]=step(Gclosed);%ess=0.29%

