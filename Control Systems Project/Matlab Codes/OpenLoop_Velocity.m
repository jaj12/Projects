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
numv=[K*n 0]
den=[L*Jo,L*bo+R*Jo,R*bo+K^2,0];
Gv=tf(numv,den)
step(Gv)
[y,t]=step(Gv);
stepinfo(Gv)