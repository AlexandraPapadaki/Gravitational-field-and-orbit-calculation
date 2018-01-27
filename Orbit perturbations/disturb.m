function [dw, dW, dM, n] = disturb(a,e,i)

a_semi = 6378.14*1000; %km
GM = 398600.44; %(km^3)(s^(-2))
J2 = 1.08263*10^(-3);

n = sqrt(GM/(a*10^(-3))^3); %h monada pou prokuptei einai s

const = 3*n*J2*a_semi^2;

dw = -const*(1 - 5*(cos(i))^2)/(4*(1-e^2)^2*a^2);

dW = -const*cos(i)/(2*(1-e^2)^2*a^2);

dM = n + const*(3*((cos(i))^2)-1)/(4*(1-e^2)^(3/2)*a^2);
end