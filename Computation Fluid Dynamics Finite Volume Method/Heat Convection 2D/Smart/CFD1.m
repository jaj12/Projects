tic
%clear; clc;
%% Solver
% Solver='GS' for Gauss Sigel
% Solver='TDMA' for TDMA
lambda=0.3;%underrelaxation coefficient
Solver="GS" ;
ItMAX= 1e6;
tol=1e-9;
ConvectionScheme= "Smart";
%Initializing P and Q for TDMA
P=0;
Q=0;

%%Intialization of geometry

L=1;%range for x-coordinates
W=L;%range for y-coordinates 

Xp=5;%number of x surfaces
Yp=5;%number of y surfaces

x=linspace(0,L,Xp); %x-geometry (in this case uniform)
y=(linspace(0,W,Yp))'; %y-geometry (in this case uniform)
%x=[0,0.1,0.3,0.55,0.6,0.8,0.9,1];
%y=[0,0.1,0.3,0.35,0.4,0.6];

[dx,dy,dxc,dyc,xc,yc,V,n,m]=geoInt(x,y); %Calculating all geometry needed variables

[gn,gs,ge,gw]= Geo(V);%geometric factors
%% input parameters
Qg=0; %heat generation W/m3

Tc=0.5*ones(n-1,m-1); %temperature matrix in Kelvin Initial Guess
Tcold=Tc;
tolMatrix=tol*ones(n-1,m-1);% Tolerence Matrix

Rho=1*ones(n-1,m-1); %Density kg/m3
Cp=1; %specific heat J/Kg. K
% "Upwind" for upwind scheme
% "Smart" for Smart scheme

%Speed Matrix
Ux=1*cosd(45)*ones(n-1,m-1);
Uy=1*cosd(45)*ones(n-1,m-1);
%% Boundry Conditions
% 1 for Constant T
% 2 for Constant Q
% 3 for Mixed
% 4 for Inlet
% 5 for Outlet
BN= 5;
BS= 4;
BE= 5;
BW= 4;

%constant Tb 
TbN= 340*ones(n-1,1);
TbS= 340*ones(n-1,1);
TbE= 340*ones(m-1,1);
TbW= 350*ones(m-1,1);


%constant Qb 
QbN= 0;
QbS= 0;
QbE= 0;
QbW= 0;

%Mixed
hbN= 0;
hbS= 0;
hbE= 0;
hbW= 0;
TinfN= 0;
TinfS= 0;
TinfE= 0;
TinfW= 0;

%Inlet
TNin=0;
TSin=0;
TEin=0;
TWin=1;


%Outlet


%% geometric loops on control volumes

%Initializing Coefficient Matrices
a=0;
ae=0;
aw=0;
an=0;
as=0;
b=0;

for k=1:ItMAX
            kc=cond(Tc);%conductivity at centers
            
            [kn,ks,ke,kw]=condsurface(kc,gn,gs,ge,gw);%conductivity at surfaces
            %calculating Border Temperatures
            [TbN,TbS,TbE,TbW]=borderTemp(Tc,TbN,TbS,TbE,TbW,BW,BE,BN,BS,QbW,QbE,QbS,QbN,hbW,hbE,hbS,hbN,TinfW,TinfE,TinfS,TinfN,dx,dy,dxc,dyc,TWin,TEin,TNin,TSin);
            
            %obtaining coefficients for the energy balance at each CV.
            [a,an,as,aw,ae,b]=CoeffCal(Tc,dxc,dyc,dy,dx,Qg,BE,BW,BN,BS,hbE,hbW,hbN,hbS,TinfE,TinfW,TinfN,TinfS,TbE,TbW,TbN,TbS,QbE,QbW,QbN,QbS,kn,ks,ke,kw,Ux,Uy,ConvectionScheme,Rho,TWin,TNin,TSin,TEin,Cp);
            
            %under relaxation
            b=b+((1-lambda)/lambda).*a.*Tcold;
            a=a/lambda;
            %Iterative Solving for Temperature at each center of CVs.
            if Solver=="GS" %Gauss-Siegel
                    [Tc]=GaussSiegel(Tc,a,an,as,ae,aw,b,n,m);
            elseif Solver== "TDMA" %TDMA
                    [Tc]=TDMA(Tc,a,an,as,ae,aw,b,n,m);
            end

    if abs(Tc-Tcold)<=tolMatrix
        it=k;
        break
    end
    %Saving temperature values for comparison.
    Tcold=Tc;
    
end


Tcomplete=[TbW,Tc',TbE'];
TbE=flip(TbE);
TbS=flip(TbS);
TbNWcorner=0.5*(TbN(1)+TbW(end));
TbNEcorner=0.5*(TbN(end)+TbE(1));
TbSWcorner=0.5*(TbS(1)+TbW(1));
TbSEcorner=0.5*(TbS(end)+TbE(end));
TNcomplete= [TbNWcorner,TbN',TbNEcorner];
TScomplete= [TbSWcorner,TbS',TbSEcorner];
Tcomplete=[TScomplete;Tcomplete;TNcomplete];

figure(2)

xplot=[0,xc,L];
yplot=[0,yc,W];
contourf(xplot,yplot,Tcomplete)
 title('Temperature(K)')
 xlabel('Length')
 ylabel('Width')
 color = colorbar; colormap(jet);
 ylabel(color, 'Temperature(K)')

%plot(yc,Tc(15,:))


toc



