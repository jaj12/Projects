  % MECH 510 project
clc; clear;

%% Initial values
%% constants
Dp = 12.7e-3; H = 0.75; Lt = 2.44;% meters
N = 100000;
syms x % variable to solve for turbulent friction factor
P = pi*Dp; Acs = pi/4*Dp^2;
Kcopper = 401; %W/m.K

Tamb = 25 + 273.15; % kelvin

%% air constants at 30 celcius
rhoAir = 1.1614; CpAir = 1.007e3; viscAir = 15.89e-6; 
KAir = 26.3e-3; alpha = 22.5e-6; PrAir = 0.707;

Th = zeros(1,10000); Tc = zeros(1,10000);

%% design variables
% Ts = 110 + 273.15; Nf = 15;
for Nf = 1:15
     Xh = linspace(0,H,10000); Xc = linspace(0,Nf*50e-3,10000);
     Af = Nf*2*30*50e-6*6; Ao = Af + Nf*(P - 6e-3)*50e-3;
    for Ts = (30:100)+273.15
        %% initial guesses
        massFlow = 2e-3; %kg/s
        Th(1) = 40+273.15;

        %% initialize the properties at T = Th(1)
        Pr = -0.0568*Th(1) +22.351; 
        % Beta = 7e-6*Ts-0.0017;
        visc = -8e-6*Th(1) + 0.0032;
        Cp= ((9e-6)*Th(1)^2 - 0.0055*Th(1) + 5.0226)*1e3; 
        Kw =  0.0009*Th(1) + 0.3353; 
        % Rho = -0.5172*Ts + 1153.6;
        Rho = -0.003*Th(1)^2 + 1.4574*Th(1) + 824.96;
        BetaAir = 1/((Tamb+Th(1))/2);
        TcAvg = Th(1);
        

        %% set old flow rate to inf to enter while loop
        massFlowOld = inf; ThinOld = inf; n = 0;
        while(abs(massFlow-massFlowOld)>1e-7 && abs(Th(1)-ThinOld) > 1e-3 && n<N)
            ThinOld = Th(1);
            massFlowOld = massFlow;
            % internal convection
            Re = 4*massFlow/(P*visc);

            if Re <= 2300 %Laminar
               Nuh = 3.66;
               Nuc = 3.66;
               f = 64/Re;
            else
                Nuh = 0.0243*Re^(4/5)*Pr^0.4;
                Nuc = 0.0265*Re^(4/5)*Pr^0.3;
                f = sym2poly(vpasolve(1/sqrt(x) == -2*log10(0.0006/3.7+2.51/(Re*sqrt(x))),x));
        %         A = (2.457*log(1/((7/Re)^0.9+0.27*0.0006)))^16;
        %         B = (37530/Re)^16;
        %         f = 8*((8/Re)^12+1/(A+B)^1.5)^(1/12);

            end
            hH = Nuh*Kw/Dp;
            hCi = Nuc*Kw/Dp;

            %natural convection
            Gr = 9.81*BetaAir*abs(TcAvg-Tamb)*(Nf*50e-3)^3/(viscAir^2);
            Ra = Gr*PrAir;
            if Ra<=1e9
                Nuco = 0.59*Ra^(1/4);
            else 
                Nuco = 0.1*Ra^(1/3);
            end
        %     Nuco = (0.825 + 0.387*Ra^(1/6)/(1+(0.492/PrAir)^(9/16))^(8/27))^2;

            hCo = Nuco*KAir/(Nf*50e-3);

            Th = Ts - (Ts - Th(1))*exp(-hH*P.*Xh/(massFlow*Cp));
        %     ThAvg = Ts + (Ts - Th(1))*exp(-hH*P*H/(massFlow*Cp))*massFlow*Cp/(P*hH*H) - (Ts - Th(1))*massFlow*Cp/(P*hH*H);
            ThAvg = trapz(Xh,Th)/H;

            m = sqrt(2*hCo/(Kcopper*1e-3));
            eta_f = tanh(m*30e-3)/(m*30e-3);
            eta_o = 1 - (1 - eta_f)*Af/Ao;
            UAf = (1/(hCi*P*Nf*50e-3)+1/(eta_o*hCo*(Nf*50e-3*(P-6e-3)+Nf*2*6*30*50e-6))).^-1;
            Uf = UAf/Af;
            Tc = Tamb - (Tamb - Th(end))*exp(-Uf*2*30e-3.*Xc*6./(massFlow*Cp));
            TcAvg = trapz(Xc,Tc)/(Nf*50e-3);
            if Th(1)-273.15>100,break,end
            Th(1) = Tc(end);

            To = (TcAvg+ThAvg)/2;

            Pr = -0.0568*To +22.351;
            Beta = 7e-6*To-0.0017;
            visc = -8e-6*To + 0.0032;
            Cp= ((9e-6)*To^2 - 0.0055*To + 5.0226)*1e3;
            Kw =  0.0009*To + 0.3353; 
        %     Rho = -0.5172*To + 1153.6;
            Rho = -0.003*To^2 + 1.4574*To + 824.96;
            BetaAir = 1/((Tamb+TcAvg)/2);

        %     RhoC = Rho*(1+Beta*(To-TcAvg));
        %     RhoH = Rho*(1+Beta*(To-ThAvg));
        %     RhoC = -0.5172*TcAvg + 1153.6;
        %     RhoH = -0.5172*ThAvg + 1153.6;
            RhoC = -0.003*TcAvg^2 + 1.4574*TcAvg + 824.96;
            RhoH = -0.003*ThAvg^2 + 1.4574*ThAvg + 824.96;
            massFlow = sqrt(abs(RhoC-RhoH)*9.81*H*2*Rho*Acs^2*Dp/(Lt*f));

            n = n + 1;
        end
        if max(Th)-273.15<100
            SteadyFlow(Ts-29-273.15,Nf) = massFlow;
            SteadyTcAvg(Ts-29-273.15,Nf) = TcAvg;
            SteadyThAvg(Ts-29-273.15,Nf) = ThAvg;
            SteadyhH(Ts-29-273.15,Nf) = hH;
        else
            SteadyFlow(Ts-29-273.15:end,Nf) = NaN;
            SteadyTcAvg(Ts-29-273.15:end,Nf) = NaN;
            SteadyThAvg(Ts-29-273.15,Nf) = NaN;
            SteadyhH(Ts-29-273.15,Nf) = NaN;
            break
        end
    end
end
