%% before control

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
G=tf(num,den)
Gv=tf(numv,den)
%Gclosed=feedback(G,1)
%step(Gclosed)
step(G)
step(Gv)
[A,B,C,D]=tf2ss(num,den)
pzmap(G)% we have a pole on the zero hence the system is marginally stable
%when we subject the system to a step input we can see that it diverges
%the system is type 1 with 3rd order
%since it has a constant slope this makes sense, since the velocity is
%constant.
step(15*G)
stepinfo(15*G)
rlocus(G)
Po=4;
Ts=2;
[zeta,wn]=SecondOrderResponse(Po,Ts)
%zeta=0.75;
zeta=0.8;
%wn=5;
sgrid(zeta,wn)
%we can't use simple P controller since this is not on rootlocus
%we want transient response of Po=4% and Ts=2s
%so, we will need PD controller
%since we also want a steady state error less than 2%
%we use PID;
% we can also use Lead instead of PD
% and lag instead of I
% so we use Lead-lag compensator

%% Pd controller

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
rlocus(PD*G)
sgrid(zeta,wn)
[Kpd poles]=rlocfind(PD*G)%Kpd=127.0620
Gpd=feedback(Kpd*PD*G,1);
step(Gpd)
stepinfo(Gpd)
[y,t]=step(Gpd);
%overshoot is 4.77% not okayso
%settling time 1.9605 which is okay

%% PID

zpid1=-5;
Dnew=D+angle(Sp(1))-angle(Sp(1)-zpid1)
zpi=-(imag(Sp(1))/tan(Dnew))+real(Sp(1))
PID=tf(conv([1 -zpid1],[1 -zpi]),[1 0])
rlocus(PID*G)
sgrid(zeta,wn)
[KpiD poles]=rlocfind(PID*G)%Kpid=169.4371
Gc2=feedback(KpiD*PID*G,1);
[y,t]=step(Gc2);
stepinfo(Gc2)
%overshoot=20.9896%
%settling time 1.7589;
%% PI

Dnew=D+angle(Sp(1))
zpi=-2;
PI=tf([1 -zpi],[1 0])
rlocus(PI*G)
sgrid(zeta,wn)
Kpi=500;
Gc2=feedback(Kpi*PI*G,1);
[y,t]=step(Gc2);
stepinfo(Gc2)
%overshoot=20.9896%
%settling time 1.7589;

%%
% Lead Compensator
zL=-P(3);
PoleAngle=atan2(imag(Sp(1)+zL),real(Sp(1)+zL))-D;
% graphically, find angle between pole and SP1, then add abs value of the
% real part of SP1 to get the pole
pL=(imag(Sp(1))/tan(PoleAngle))-real(Sp(1));
Gc=tf([1 zL],[1 pL]);
T1=1/zL;
gamma=pL*T1;
rlocus(Gc*G)
sgrid(zeta,wn)
% axis([-6 1 -10 10]);
[Kc,pol]=rlocfind(Gc*G);%750.3234 for z=threshold not enough settling time is high,K=746.808 for zeta=0.75
Gc3=feedback(Kc*Gc*G,1);
step(Gc3)
stepinfo(Gc3)%K=751.0046 for z=0.78 w=diff we are happy with these., 
[y,t]=step(Gc3);
%% Lag Compensator
% We know that we use the lag compensator to fix our steady state response,
% when we are satisfied with our transient response, since this is not the
% case for us, we will not even consider a stand-alone Lag compensator.
%% Lead-Lag Compensator
% we already have an acceptable steady state error for our system in the
% lead compensator, so we do not need to add a lag compensator with it.
num3=2.0959;
GL=tf(num3,den)
bode(GL)
margin(GL)
