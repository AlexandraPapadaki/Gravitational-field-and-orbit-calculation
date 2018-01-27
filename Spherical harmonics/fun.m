function [Cnm,Snm]=fun(nn,mm,theta,lambda)
         Pnm = legendre (nn,mm,theta);
         Cnm = Pnm*cos(mm*lambda*pi/180);
         Snm = Pnm *sin(mm*lambda*pi/180);
end

