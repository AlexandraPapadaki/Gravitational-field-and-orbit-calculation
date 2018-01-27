clear all
close all
clc


%data input
GM=398600.44;
w_earth=2*pi()/86164;
data=load('data.txt');
data(1,3:5)= data(1,3:5)*pi()/180;
M = 64.0942*pi/180;

Eo = M;
E(1)=M+data(1,2)*sin(Eo);
for k=1:50
    E(k+1)=M+data(1,2)*sin(E(k));
    if abs(E(k+1)-E(k))<10e-8
        break
    end
end

Ec_rad = E(:,end);

K=[cos(data(1,4));sin(data(1,4));0];
C0=[sin(data(1,4))*sin(data(1,3));-cos(data(1,4))*sin(data(1,3));cos(data(1,3))];
P=cos(data(1,5))*K(:,1)+sin(data(1,5))*(cross(C0(:,1),K(:,1)));
Q=-sin(data(1,5))*K(:,1)+cos(data(1,5))*(cross(C0(:,1),K(:,1)));

f_rad=2*atan((tan(Ec_rad/2)*sqrt((1+data(1,2))/(1-data(1,2)))));
r=data(1,1)*(1-data(1,2)*cos(Ec_rad));
p=r*(1+data(1,2)*cos(f_rad));
c=sqrt(p*GM);

r_pos(1,1:3)=r*(cos(f_rad)*P+sin(f_rad)*Q);

vel(1,1:3)=(c/p)*(-sin(f_rad)*P+(data(1,2)+cos(f_rad))*Q);
    
x_initial = [vel r_pos];

R_earth = 6378.137;
T_orbit=2*pi*sqrt((data(1,1)^3)/GM);

%% Numerical itegration
v=sqrt(GM/data(1,1));
tspan =(0:1:T_orbit*3)';
y0 = x_initial;
% y0(3,1) = v;
% y0(4,1) = data(1,1);

options12 = odeset('RelTol',1e-12,'AbsTol',[1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]);
[T,ys113]=ode113(@yprime_per,tspan,y0',options12);

for i=1:length(ys113)
  r_new(i,1) = sqrt(ys113(i,4)^2+ys113(i,5)^2+ys113(i,6)^2);  
  h(i,1) = r_new(i,1) - R_earth;
  
  % the magnitude of the velocity in the Earth-fixed system
  B = [0; 0; w_earth];
  x_atm(1:3,i) = cross(B,ys113(i,4:6)');
  dif_x(i,1:3) = ys113(i,1:3)-x_atm(1:3,i)';
  velocity(i,1) = sqrt(dif_x(i,1)^2+dif_x(i,2)^2+dif_x(i,3)^2);
end

%Calculation corrections because of earths rotation
E_int=1*360/T_orbit;

pos_true = zeros(length(ys113),3);
for i=1:length(ys113)
    
    t_dE=(i)*T_orbit*E_int/360;
    dw_ang=w_earth*t_dE;
    R_3=[cos(-dw_ang) -sin(-dw_ang) 0;sin(-dw_ang) cos(-dw_ang) 0;0 0 1];
    pos_true(i,1:3)= ys113(i,4:6)*R_3';

end

%% Plots
plot(1:15141,velocity(1:15141,1));
xlabel('time(sec)')
ylabel('Magnitude of the velocity in the Earth-fixed system')
grid on 
figure
plot(1:15141,h(1:15141,1));
xlabel('time(sec)')
ylabel('Height of the satellite above ground')
grid on 

%%
figure
earth_sphere
hold on
x = pos_true(1:15141,1);
y = pos_true(1:15141,2);
z = pos_true(1:15141,3);
plot3(x,y,z,'.r');