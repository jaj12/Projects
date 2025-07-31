tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note that the code didnot run. However, all the needed steps and equations
%were added. in some cases I did not take the boundary conditions into
%account like in the RhieChow interpolation for example and Conituity
%pressure corection equation. 

%The steps left are: 
% 1) fixing the Rhiechow interpolation and Pressure
%correction equation such that boundary control volumes are taken into
%account.

% 2) fixing the needed input arguements for the VelocityX and VelocityY
%function.

%Adding cases 5 and 6 to the boundary conditions.

%I honestly believe if I started my code off since assignment 1 knowing I
%will need to incorporate velocity boundary conditions like this it would
%have been a cleaner code and easier to integrate to this assignment.

%I know the "Temperature" function should've been sufficient to use for
%VelocityX and VelocityY Initially, but unfortunetly the boundary
%conditions had to all change.

%Thank you for a very enjoyable course and learning experience!.





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear; clc;
%% Solver
% Solver='GS' for Gauss Sigel
% Solver='TDMA' for TDMA
lambda=1;%underrelaxation coefficient
%underrelaxation coefficients for temperature, velocity and pressure
%correction should be done seperately such the the sum of pressure lambda
%and velocity lambda is equal to 1.1
Solver="GS" ;
ItMAX= 1e5;%iteration max for convergence of solution
tol=1e-6;%tolerance for convergence of solution (Residual)
Ptol=1e-6; %Pressure tolerence for convergence
Pitmax=10000; %iteration max for pressure correction and navier stokes
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
TransientScheme="none";

%Initializing P and Q for TDMA
P=0;
Q=0;

%Intialization of geometry

L=2;%range for x-coordinates
W=1;%range for y-coordinates 

Xp=41;%number of x surfaces
Yp=21;%number of y surfaces

x=linspace(0,L,Xp); %x-geometry (in this case uniform)
y=(linspace(0,W,Yp))'; %y-geometry (in this case uniform)
%x=[0,0.1,0.3,0.55,0.6,0.8,0.9,1];
%y=[0,0.1,0.3,0.35,0.4,0.6];

[dx,dy,dxc,dyc,xc,yc,V,n,m]=geoInt(x,y); %Calculating all geometry needed variables

[gn,gs,ge,gw]= Geo(V);%geometric factors



%% input parameters
Qg=0; %heat generation W/m3
Fb=0; %body forces N
%Temperature Matrix Initial Guess
Tc=1000*ones(n-1,m-1); %Kelvin
Tcold=Tc;
Tcoldtrans=Tc;
%Speed Matrix Initial Guess
Ux=0*ones(n-1,m-1); %m/s
Uy=0*ones(n-1,m-1); %m/s
Uxoldtrans=Ux;
Uyoldtrans=Uy;
%Pressure Intial Guess
P=10000*ones(n-1,m-1);
Pcorr=10e-6*ones(n-1,m-1);%Pressure correction initialiation.
%Pressure at boundaries
PbN=10000*ones(n-1,1);
PbS=10000*ones(n-1,1);
PbE=10000*ones(1,m-1);
PbW=10000*ones(1,m-1);
%
Rho=0.8*ones(n-1,m-1); %Density kg/m3
Cp=1030; %specific heat J/Kg. K
%conductivity is in cond.m
%viscosity is in visco.m

%% Boundry Conditions
% 1 for Constant T
% 2 for Constant Q
% 3 for Mixed
% 4 for Inlet Specified Velocity
% 5 for Inlet Specified Static Pressure and Velocity Direction
% 6 for Inlet Specified Total Pressure and Velocity Direction
% 7 for Outlet Specified Static Pressure
% 8 for Outlet Fully Developed Outlet Flow
% 9 for Wall Boundry No Slip
% 10 for Wall Boundry Slip


%the boundaries are divided into 2 sections 1 and 2;
% first Boundary condition is for the 1st section 
% and the 2nd is for the second section

% North
WN=[9 0];% Wall condition

BN= [1 4];
N1=1.8;
N2=L-N1;
n2i=n-ceil(N2/L*(n-1));
% South
WS=[9 9];

BS= [1 1];
S1=2;
S2=L-S1;
s2i=n-ceil(S2/L*(n-1));
% East
WE=[9 9];

BE= [2 2];
E1=1;
E2=W-E1;
e2i=m-ceil(E2/W*(m-1));
% West
WW=[0 9];

BW= [7 3];
W1=0.1;
W2=W-W1;
w2i=m-ceil(W2/W*(m-1));
% 1 constant T
TbN= 300*ones(n-1,1);
TbS= 500*ones(n-1,1);
TbE= 0*ones(m-1,1);
TbW= 0*ones(m-1,1);


% 2 constant Q
QbN= 0;
QbS= 0;
QbE= 0;
QbW= 0;

% 3 Mixed
hbN= 0;
hbS= 0;
hbE= 0;
hbW= 15;
TinfN= 0;
TinfS= 0;
TinfE= 293;
TinfW= 0;

% Inlet 
TNin=0;
TSin=0;
TEin=0;
TWin=1;

% 4 Specified Velocity
UEinx=0;
UWinx=0;
UNinx=5;
USinx=0;

UEiny=0;
UWiny=0;
UNiny=5;
USiny=0;
% 5 or 6

UdirEin=0;
UdirWin=0;
UdirNin=0;
UdirSin=0;

% 5 Inlet Specified Static Pressure and Velocity Direction
PNin=0;
PSin=0;
PEin=0;
PWin=0;

% 6 Inlet Specified Total Pressure and Velocity Direction
PtotN=0;
PtotS=0;
PtotE=0;
PtotW=0;

% Outlet

% 7 Outlet Specified Static Pressure

PNout=0;
PSout=0;
PEout=0;
PWout=1e5;

% 8 Outlet Fully Developed Outlet Flow

% 9 for Wall Boundry No Slip
UWwallx=0;
UNwallx=0;
USwallx=0;
UEwallx=0;

UWwally=0;
UNwally=0;
USwally=0;
UEwally=0;
% 10 for Wall Boundry Slip

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
                muc=visco(Tc);%viscosity at centers
                
                for i=1:Pitmax
                [UbWx,UbWy,UbEx,UbEy,UbSx,UbSy,UbNx,UbNy]= boundaryvelocity(UWwallx,UWwally,UEwallx,UEwally,UNwallx,UNwally,USwallx,USwally,BE,BW,BN,BS,WN,WE,WS,WW,Ux,Uy,UNiny,USiny,UNinx,USinx,UEinx,UWinx,n2i,s2i,e2i,w2i,n,m);
                [Uwx,Uex,Uny,Usy]=speedsurface(Ux,Uy,gn,gs,ge,gw);%Surface Velocities
                [kn,ks,ke,kw]=condsurface(kc,gn,gs,ge,gw);%conductivity at surfaces
                [mun,mus,mue,muw]=viscosurface(muc,gn,gs,ge,gw);%viscosity at surfaces
                [Pn,Ps,Pe,Pw]=Pressuresurface(P,gn,gs,ge,gw);%Pressure at surfaces
                %obtaining coefficients for the momentum equation.
                %some input arguements to the VelocityX and VelocityY
                %functions might be missing.
                [awx,anx,aex,asx,ax,bx,dpcxx,dpcyx]=VelocityX(P,dxc,dyc,dy,dx,BE,BW,BN,BS,WN,WE,WS,WW,Ux,Uy,ConvectionScheme,Rho,dt,TransientScheme,UNiny,USiny,UNinx,USinx,UEinx,UWinx,Uwx,Uex,Uny,Usy,Uxoldtrans,TbN,TbS,TbE,TbW,n2i,s2i,e2i,w2i,PbW,PbE,PbS,PbN);
                [awy,any,aey,asy,ay,by,dpcxy,dpcyy]=VelocityY(P,dxc,dyc,dy,dx,BE,BW,BN,BS,WN,WE,WS,WW,Ux,Uy,ConvectionScheme,Rho,dt,TransientScheme,UNiny,USiny,UWiny,UEiny,UEinx,UWinx,Uwx,Uex,Uny,Usy,Uyoldtrans,TbN,TbS,TbE,TbW,n2i,s2i,e2i,w2i,PbW,PbE,PbS,PbN);
                %Solving for Ux and Uy
                
                %underrelaxation
                bx=bx+((1-lambda)/lambda).*ax.*Uxold;
                ax=ax/lambda;
                
                by=by+((1-lambda)/lambda).*ay.*Uyold;
                ay=ay/lambda;
                if Solver=="GS" %Gauss-Siegel
                        [Ux]=GaussSiegel(Ux,ax,anx,asx,aex,awx,bx,n,m);
                        [Uy]=GaussSiegel(Uy,ay,any,asy,aey,awy,by,n,m);
                elseif Solver== "TDMA" %TDMA
                        [Ux]=TDMA(Ux,ax,anx,asx,aex,awx,bx,n,m);
                        [Uy]=TDMA(Uy,ay,any,asy,aey,awy,by,n,m);
                end
                %RhieChow Interpolation to find mdot*
                [mdotn,mdotw,mdote,mdots]=RhieChow(ac,P,dxc,dyc,dx,dy,ge,gn,Rho,dpcxx,dpcyy,Ux,Uy);
                %Pressure correction
                [acp,anp,asp,aep,awp,bp]=Continuity(dxc,dyc,dy,dx,Ux,Uy,Rho,mdotw,mdote,mdotn,mdots,acx,acy,n,m);
                %Solving for Pcorr
                if Solver=="GS" %Gauss-Siegel
                        [Pcorr]=GaussSiegel(Pcorr,acp,anp,asp,aep,awp,bp,n,m);                  
                elseif Solver== "TDMA" %TDMA
                        [Pcorr]=TDMA(Up,ap,anp,asp,aep,awp,bp,n,m);
                end
                if Max(Pcorr)<=Ptol
                    break
                end
                P=P+Pcorr;%This might need explicit under relaxation to avoid divergence of pressure.
                end
                %calculating Border Temperatures
                [TbN,TbS,TbE,TbW]=borderTemp(Tc,TbN,TbS,TbE,TbW,BW,BE,BN,BS,QbW,QbE,QbS,QbN,hbW,hbE,hbS,hbN,TinfW,TinfE,TinfS,TinfN,dx,dy,dxc,dyc,TWin,TEin,TNin,TSin);

                %obtaining coefficients for the energy balance.
                [a,an,as,aw,ae,b]=Temperature(Tc,dxc,dyc,dy,dx,Qg,BE,BW,BN,BS,hbE,hbW,hbN,hbS,TinfE,TinfW,TinfN,TinfS,TbE,TbW,TbN,TbS,QbE,QbW,QbN,QbS,kn,ks,ke,kw,Ux,Uy,ConvectionScheme,Rho,TWin,TNin,TSin,TEin,Cp,dt,TransientScheme,Tcoldtrans,UNyin,USyin,UExin,UWxin,Uwx,Uex,Uny,Usy);
                
              
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
        %Saving previous velocity values.
        Uyold=Uy;
        Uxold=Ux;
    end
            if TransientScheme=="CN"
            Tc=2.*Tc-Tcoldtrans;
            end
    
    if max(abs(Tc-Tcoldtrans))<=toltrans
            it=k;
            break
    end
   
    Tcoldtrans=Tc;
    Uxoldtrans=Ux;
    Uyoldtrans=Uy;
end




figure(1)
contourf(xc,yc,Tc')
title('Temperature(K)')
xlabel('Length')
ylabel('Width')
color = colorbar; colormap(jet);
ylabel(color, 'Temperature(K)')

figure(2)
contourf(xc,yc,P')
title('Pressure (Pa)')
xlabel('Length')
ylabel('Width')
color = colorbar; colormap(jet);
ylabel(color,'Pressure (Pa)')

figure(3)
contourf(xc,yc,Ux')
title('Velocity in x (m/s)')
xlabel('Length')
ylabel('Width')
color = colorbar; colormap(jet);
ylabel(color, 'Velocity in x (m/s)')

figure(4)
contourf(xc,yc,Uy')
title('Velocity in y (m/s)')
xlabel('Length')
ylabel('Width')
color = colorbar; colormap(jet);
ylabel(color, 'Velocity in y (m/s)')




toc



