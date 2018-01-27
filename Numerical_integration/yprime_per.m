function ys = yprime_per(t,y0)
          GM = 398600.44;
          CD = 3;
          A = 0.954/10^6;
          m = 872;
          w_earth=2*pi()/86164;
          R_earth = 6378.137;
          ys = zeros(6,1);
                   
          %velocity of the atmosphere 
          B = [0; 0; w_earth];
          x_atm = cross(B,y0(4:6));
          
          % difference between the velocity of the satellite and velocity
          % of the atmosphere
          dif_x = y0(1:3)-x_atm;
          r_dif = sqrt(dif_x(1)^2+dif_x(2)^2+dif_x(3)^2);
                            
          %density of the atmosphere
          r = sqrt(y0(4)^2+y0(5)^2+y0(6)^2);
          h = r - R_earth;
          ro = 1.025^(500-h)/1000;
          
          % atmospheric perturbations
          x_D = -1/2*CD*ro*A/m*r_dif*dif_x;
          
          ys(1) = -GM/(r^3)*y0(4)+x_D(1);
          ys(2) = -GM/(r^3)*y0(5)+x_D(2);
          ys(3) = -GM/(r^3)*y0(6)+x_D(3);
          ys(4) = y0(1);
          ys(5) = y0(2);
          ys(6) = y0(3);
end
