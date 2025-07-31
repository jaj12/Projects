close all; clear; clc;
Jm=0.04;
Jl=0;
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
Kv=2.43;
K=629.39
G1=tf([2.0959],den)
Pneeded=24.2;
alpha=0.4185;
w=2.67;
T=0.578;
z=1.727;
p=4.12763;
Kc=1503.918757;
Gc=tf([1 z],[1 p]);
Gfull=feedback(Kc*Gc*G,1);
step(Gfull)
stepinfo(Gfull)
margin(Kc*Gc*G)


