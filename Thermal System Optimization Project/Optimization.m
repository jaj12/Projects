clc; clear;
% Mech 510 Optimization

%% constants
Dp = 12.7e-3; H = 0.75; Lt = 2.44;% meters
N = 100000;
syms x % variable to solve for turbulent friction factor
P = pi*Dp; Acs = pi/4*Dp^2;
Kcopper = 401; %W/m.K

Tamb = 25 + 273.15; % kelvin
%% objective function and constraints
interest = 0.0067;
interestVector = 1./(1+interest).^(1:12);
interestSum = sum(interestVector);

%X(1) = Ts, X(2) = Nf 

massFlow = @(X) (-0.00086618 + -4.0304e-06*X(1) + 2.4422e-08*X(1)^2) + (-6.6765e-05 + -4.9464e-07*X(1) + 2.8121e-09*X(1)^2)*X(2) + (-3.6277e-06 + 4.6695e-08*X(1) + -1.3234e-10*X(1)^2)*X(2)^2;
TcAvg    = @(X) (2.5974 + 0.99057*X(1) + 4.6845e-07*X(1)^2) + (6.021 + -0.024819*X(1) + 1.564e-05*X(1)^2)*X(2) + (-0.29418 + 0.0015591*X(1) + -1.9251e-06*X(1)^2)*X(2)^2;
ThAvg    = @(X) (4.3383 + 0.98404*X(1) + 6.4927e-06*X(1)^2) + (4.5568 + -0.016783*X(1) + 4.9266e-06*X(1)^2)*X(2) + (-0.22922 + 0.0012198*X(1) + -1.4939e-06*X(1)^2)*X(2)^2;
hH       = @(X) (102.7897 + 0.22485*X(1) + 4.6706e-05*X(1)^2) + (-0.28775 + 0.004837*X(1) + -1.3045e-05*X(1)^2)*X(2) + (0.048527 + -0.00035445*X(1) + 6.4776e-07*X(1)^2)*X(2)^2;

objectiveFunction = @(X) (1.35*X(2) + 48*(10^-3)*hH(X)*P*H*(X(1)-ThAvg(X))*interestSum);
%%
%constraints
G = @(X)([1.5 - 10^3*massFlow(X), -308.15 + TcAvg(X)]);
H = @(X) (0);
s = 1; w = 1;
LB = [303.15, 1]; UB = [373.15, 15];
X0 = [350, 4];
tolh = 1e-3; tolg = 1e-4; tolCyclic = 1e-8; tolGolden = 1e-8;
[Xopt, i] = Penalty(objectiveFunction, G, H, s, w, LB, UB, X0, tolh, tolg, tolCyclic, tolGolden);

Xopt(2) = ceil(Xopt(2));
f = @(Ts)(objectiveFunction([Ts,Xopt(2)]));
G = @(Ts)(G([Ts,Xopt(2)]));
H = @(Ts) (0);
s = 1; w = 1;
Xopt(1) = Penalty(f,G, H, s, w, 303.15, 373.15, 350, tolh, tolg, tolCyclic, tolGolden);
Xopt
