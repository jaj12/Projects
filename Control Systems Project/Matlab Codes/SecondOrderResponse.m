function [zeta wn] = SecondOrderResponse(PO,Ts)
A = (log(PO/100))^2;
zeta = sqrt(A/(A+(pi^2)));
wn = 4/(zeta*Ts);