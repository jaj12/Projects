tic
%clear; clc;
%% Solver
% Solver='GS' for Gauss Sigel
% Solver='TDMA' for TDMA
lambda=1;%underrelaxation coefficient
Solver="GS" ;
ItMAX= 1e5;%iteration max for convergence of solution
tol=1e-6;%tolerance for convergence of solution (Residual)

% "Upwind" for upwind scheme
% "Smart" for Smart scheme
ConvectionScheme= "Smart";

% Transient 
%in seconds
ItMAXtrans= 60;%max time iterations
dt=1;
toltrans=1e-6;%tolerance for reaching steady state
%"BE" for backward euler
%"CN" for crank nicolson
TransientScheme="CN";

%Initializing P and Q for TDMA
P=0;
Q=0;

%Intialization of geometry

L=0.1;%range for x-coordinates
W=0.2;%range for y-coordinates 

Xp=21;%number of x surfaces
Yp=41;%number of y surfaces

x=linspace(0,L,Xp); %x-geometry (in this case uniform)
y=(linspace(0,W,Yp))'; %y-geometry (in this case uniform)
%x=[0,0.1,0.3,0.55,0.6,0.8,0.9,1];
%y=[0,0.1,0.3,0.35,0.4,0.6];

[dx,dy,dxc,dyc,xc,yc,V,n,m]=geoInt(x,y); %Calculating all geometry needed variables

[gn,gs,ge,gw]= Geo(V);%geometric factors



%% input parameters
Qg=0; %heat generation W/m3

Tc=1000*ones(n-1,m-1); %temperature matrix in Kelvin Initial Guess
Tcold=Tc;
Tcoldtrans=Tc;

Rho=7800*ones(n-1,m-1); %Density kg/m3
Cp=500; %specific heat J/Kg. K


%Speed Matrix
Ux=0*ones(n-1,m-1);
Uy=0*ones(n-1,m-1);
%% Boundry Conditions
% 1 for Constant T
% 2 for Constant Q
% 3 for Mixed
% 4 for Inlet
% 5 for Outlet
BN= 1;
BS= 1;
BE= 2;
BW= 2;

%constant Tb 
TbN= 0*ones(n-1,1);
TbS= 0*ones(n-1,1);
TbE= 0*ones(m-1,1);
TbW= 0*ones(m-1,1);


%constant Qb 
QbN= 0;
QbS= 0;
QbE= 0;
QbW= 0;

%Mixed
hbN= 0;
hbS= 0;
hbE= 5;
hbW= 0;
TinfN= 0;
TinfS= 0;
TinfE= 293;
TinfW= 0;

%Inlet
TNin=0;
TSin=0;
TEin=0;
TWin=1;


UExin=0;
UWxin=0;

UNyin=0;
USyin=0;

%Outlet


%% geometric loops on control volumes

%Initializing Coefficient Matrices
a=0;
ae=0;
aw=0;
an=0;
as=0;
b=0;
for l=1:ItMAXtrans
    for k=1:ItMAX
                kc=cond(Tc);%conductivity at centers
                [Uwx,Uex,Uny,Usy]=speedsurface(Ux,Uy,gn,gs,ge,gw);
                [kn,ks,ke,kw]=condsurface(kc,gn,gs,ge,gw);%conductivity at surfaces
                %calculating Border Temperatures
                [TbN,TbS,TbE,TbW]=borderTemp(Tc,TbN,TbS,TbE,TbW,BW,BE,BN,BS,QbW,QbE,QbS,QbN,hbW,hbE,hbS,hbN,TinfW,TinfE,TinfS,TinfN,dx,dy,dxc,dyc,TWin,TEin,TNin,TSin);

                %obtaining coefficients for the energy balance at each CV.
                [a,an,as,aw,ae,b]=CoeffCal(Tc,dxc,dyc,dy,dx,Qg,BE,BW,BN,BS,hbE,hbW,hbN,hbS,TinfE,TinfW,TinfN,TinfS,TbE,TbW,TbN,TbS,QbE,QbW,QbN,QbS,kn,ks,ke,kw,Ux,Uy,ConvectionScheme,Rho,TWin,TNin,TSin,TEin,Cp,dt,TransientScheme,Tcoldtrans,UNyin,USyin,UExin,UWxin,Uwx,Uex,Uny,Usy);

                %under relaxation
                b=b+((1-lambda)/lambda).*a.*Tcold;
                a=a/lambda;
                %Iterative Solving for Temperature at each center of CVs.
                if Solver=="GS" %Gauss-Siegel
                        [Tc]=GaussSiegel(Tc,a,an,as,ae,aw,b,n,m);
                elseif Solver== "TDMA" %TDMA
                        [Tc]=TDMA(Tc,a,an,as,ae,aw,b,n,m);
                end
            
        %scaled residual calculation
        [Rscaled,ss]=ResidCalc(a,aw,an,as,ae,Tc,b,n,m);
        %Condition for convergence
        if Rscaled<=tol
            it=k;
            break
        end
        %Saving temperature values for comparison.
        Tcold=Tc;
    end
            if TransientScheme=="CN"
            Tc=2.*Tc-Tcoldtrans;
            end
    
    if max(abs(Tc-Tcoldtrans))<=toltrans
            it=k;
            break
    end
    if l==10
        Tc10=Tc;
    elseif l==30
        Tc30=Tc;
    elseif l==60
        Tc60=Tc;
    end
    Tcoldtrans=Tc;
end




figure(1)
contourf(xc,yc,Tc30')
title('Temperature(K) t=30s')
xlabel('Length')
ylabel('Width')
color = colorbar; colormap(jet);
ylabel(color, 'Temperature(K)')

figure(2)
contourf(xc,yc,Tc10')
title('Temperature(K) t=10s')
xlabel('Length')
ylabel('Width')
color = colorbar; colormap(jet);
ylabel(color, 'Temperature(K)')

figure(3)
contourf(xc,yc,Tc60')
title('Temperature(K)t=60s')
xlabel('Length')
ylabel('Width')
color = colorbar; colormap(jet);
ylabel(color, 'Temperature(K)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To compare with exact solution run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Texact at tplot]=Exactsol(maxt,tplot)
[Texact60]=Exactsol(60,60,4);
[Texact30]=Exactsol(30,30,5);
[Texact10]=Exactsol(10,10,6);


toc



