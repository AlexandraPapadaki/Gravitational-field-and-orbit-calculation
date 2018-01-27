clear all
close all
clc

%data input
GM=398600.44*10^9;
w=2*pi/86164;
data=load('data_ham.txt');
data(:,3:5)= data(:,3:5)*pi()/180;

%% orbital period
T=zeros(5,1);
for i=1:4
  T(i,1)=2*pi*sqrt((data(i,1)^3)/GM);
end
T(5,1)=2*pi()/w;  %%%%%%%%  
data(5,1)=nthroot(GM*86164^2/(4*pi()^2),3);

%% position and velocity
E_int=1;
% nr_T=2;
% Ec_deg=(0:E_int:360*nr_T)';
n=7200;
Ec_deg=(0:7200)';
% Ec_rad=Ec_deg*pi/180;
Ec_rad = linspace(0,4*pi,n)';

sz=size(Ec_rad);

f_rad=zeros(sz(1,1),5);
r=zeros(sz(1,1),5);
p=zeros(sz(1,1),5);
c=zeros(sz(1,1),5);
K=zeros(3,5);
C0=zeros(3,5);
P=zeros(3,5);
Q=zeros(3,5);
r_pos=zeros(sz(1,1),15);
v=zeros(sz(1,1),15);
for i=1:5
    K(:,i)=[cos(data(i,4));sin(data(i,4));0];
    C0(:,i)=[sin(data(i,4))*sin(data(i,3));-cos(data(i,4))*sin(data(i,3));cos(data(i,3))];
    P(:,i)=cos(data(i,5))*K(:,i)+sin(data(i,5))*(cross(C0(:,i),K(:,i)));
    Q(:,i)=-sin(data(i,5))*K(:,i)+cos(data(i,5))*(cross(C0(:,i),K(:,i)));
    for j=1:sz(1,1)
        f_rad(j,i)=2*atan((tan(Ec_rad(j,1)/2)*sqrt((1+data(i,2))/(1-data(i,2)))));
        r(j,i)=data(i,1)*(1-data(i,2)*cos(Ec_rad(j,1)));
        p(j,i)=r(j,i)*(1+data(i,2)*cos(f_rad(j,i)));
        c(j,i)=sqrt(p(j,i)*GM);
        
        r_pos(j,3*(i-1)+1:3*(i-1)+3)=r(j,i)*(cos(f_rad(j,i))*P(:,i)+sin(f_rad(j,i))*Q(:,i));
        
        v(j,3*(i-1)+1:3*(i-1)+3)=(c(j,i)/p(j,i))*(-sin(f_rad(j,i))*P(:,i)+(data(i,2)+cos(f_rad(j,i)))*Q(:,i));
    end
end

%Calculation corrections because of earths rotation
pos_true = zeros(sz(1,1),15);
for i=1:5
    for j=1:7200
        t_dE(j,i)=(Ec_rad(j,1)-data(i,2)*sin(Ec_rad(j,1)))*T(i,1)/(2*pi);
        dw_ang(j,i)=t_dE(j,i)*w;
        R_3=[cos(-dw_ang(j,i)) -sin(-dw_ang(j,i)) 0;sin(-dw_ang(j,i)) cos(-dw_ang(j,i)) 0;0 0 1];
        pos_true(j,3*(i-1)+1:3*(i-1)+3)=(r_pos(j,3*(i-1)+1:3*(i-1)+3)*R_3');
    end 
end

earth_sphere
hold on
%% ERSA II 
x1 = pos_true(:,1)/1000;
y1 = pos_true(:,2)/1000;
z1 = pos_true(:,3)/1000;
plot3(x1,y1,z1,'.r');
hold on
% GPS
x2 = pos_true(:,4)/1000;
y2 = pos_true(:,5)/1000;
z2 = pos_true(:,6)/1000;
plot3(x2,y2,z2,'.k');
hold on
% CHAMP
x3 = pos_true(:,7)/1000;
y3 = pos_true(:,8)/1000;
z3 = pos_true(:,9)/1000;
plot3(x3,y3,z3,'.y');
hold on
% Molniya
x4 = pos_true(:,10)/1000;
y4 = pos_true(:,11)/1000;
z4 = pos_true(:,12)/1000;
plot3(x4,y4,z4,'.m');
hold on
% Beidou
x5 = pos_true(:,13)/1000;
y5 = pos_true(:,14)/1000;
z5 = pos_true(:,15)/1000;
plot3(x5,y5,z5,'.g');
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
hold on
grid on 
legend('','ErSa II', 'GPS', 'CHAMP', 'Molniya', 'Beidou')

%% ellips
f=1/298.257222101;
a=6378137;
e = sqrt(2*f-f^2);
wgs84=[a,e];
lat=zeros(7200,5);
lon=zeros(7200,5);
h=zeros(7200,5);

for i=1:5
    for j=1:7200
%         [lat(j,i),lon(j,i),h(j,i)] = ecef2geodetic(pos_true(j,3*(i-1)+1),pos_true(j,3*(i-1)+2),pos_true(j,3*(i-1)+3),wgs84);
        [lat(j,i),lon(j,i),h(j,i)] = xyz2blh(pos_true(j,3*(i-1)+1),pos_true(j,3*(i-1)+2),pos_true(j,3*(i-1)+3),'grs80');
    end
end
% lat=lat*180/pi;
% lon=lon*180/pi;

load('coast0.mat')
figure
plot(lam,phi)
hold on
ylabel('Latitude [deg]');
xlabel('Longitude [deg]');
title('Satellite Ground Track');
xx1 = lat(:,1);
yy1 = lon(:,1);
plot(yy1,xx1,'.r');
hold on
xx2 = lat(:,2);
yy2 = lon(:,2);
plot(yy2,xx2,'.k');
hold on
xx3 = lat(:,3);
yy3 = lon(:,3);
plot(yy3,xx3,'.y');
hold on
xx4 = lat(:,4);
yy4 = lon(:,4);
plot(yy4,xx4,'.m');
hold on
xx5 = lat(:,5);
yy5 = lon(:,5);
plot(yy5,xx5,'.g','linewidth',2);
hold on
legend('','ErSa II', 'GPS', 'CHAMP', 'Molniya', 'Beidou')

%% 3D distances
dist=zeros(7200,5);
X_ber=3783.26649*1000;
Y_ber=901.6469*1000;
Z_ber=5038.24814*1000;
pos_ber = [X_ber Y_ber Z_ber];

M2 = [1  0  0 ; 0  -1  0 ; 0  0  1];
% l_local=atan2(Y_ber,X_ber);
% f_local=atan(Z_ber/sqrt(X_ber^2+Y_ber^2));
[f_local,l_local,h_local] = xyz2blh(X_ber,Y_ber,Z_ber,'grs80');
f_local = f_local*pi/180;%%%%%
l_local = l_local*pi/180;%%%%%

rotations=[-sin(l_local),cos(l_local),0;
           -sin(f_local)*cos(l_local),-sin(f_local)*sin(l_local),cos(f_local);
            cos(f_local)*cos(l_local),cos(f_local)*sin(l_local),sin(f_local)];

        
time = zeros(7200,5);
El = zeros(7200,5);
Az = zeros(7200,5);
Zen = zeros(7200,5);
for i=1:5
    for j=1:7200
%         rot_mat = R2(pi/2-f_local)*M2*R3(l_local-pi);
%         dx(j,3*(i-1)+1:3*(i-1)+3)= pos_true(j,3*(i-1)+1:3*(i-1)+3)-pos_ber;
%         R_local(j,3*(i-1)+1:3*(i-1)+3) = (rot_mat*dx(j,3*(i-1)+1:3*(i-1)+3)')';
%         dist(j,i) = norm(R_local(j,3*(i-1)+1:3*(i-1)+3)); 
%          
%         R_local_zen_Azi = zen_azi(R_local(j,3*(i-1)+1:3*(i-1)+3)'); 
        time(j,i) = j*T(i,1)/3600;
        
        dist(j,i)=sqrt((pos_true(j,3*(i-1)+1)-X_ber)^2+(pos_true(j,3*(i-1)+2)-Y_ber)^2+(pos_true(j,3*(i-1)+3)-Z_ber)^2);
        Dr = [pos_true(j,3*(i-1)+1)-X_ber,pos_true(j,3*(i-1)+2)-Y_ber,pos_true(j,3*(i-1)+3)-Z_ber]; 
        ENU=rotations*Dr'; 

%          Zenith(j,i)=R_local_zen_Azi(1);
%          Azimuth(j,i)=R_local_zen_Azi(2);

         Elevation=(asin(ENU(3)/norm(ENU)));
         Azimuth=(atan2(ENU(1)/norm(ENU),ENU(2)/norm(ENU)));
         El(j,i)=Elevation;
         Az(j,i)=Azimuth;
         Zen(j,i)=(pi/2)-El(j,i);
    end 
end

azim = Az*180/pi();
elev = El*180/pi();
zen = Zen*180/pi();
figure

skyplot(azim(:,1),elev(:,1),'.r')
title('Skyplot - ErSa II')

figure;
skyplot(azim(:,2),elev(:,2),'.k')
title('Skyplot - GPS')

figure;
skyplot(azim(:,3),elev(:,3),'.y')
title('Skyplot - CHAMP')

figure;
skyplot(azim(:,4),elev(:,4),'.m')
title('Skyplot - Molniya')

figure;
skyplot(azim(:,5),elev(:,5),'.g')
title('Skyplot - Beidou')

figure;
ErSaII =[azim(:,1), elev(:,1), Zen(:,1), dist(:,1), time(:,1)];
GPS =[azim(:,2), elev(:,2), Zen(:,2), dist(:,2), time(:,2)];
CHAMP = [azim(:,3), elev(:,3), Zen(:,3), dist(:,3), time(:,3)];
Molniya =[azim(:,4), elev(:,4), Zen(:,4), dist(:,4), time(:,4)];
Beidou =[azim(:,5), elev(:,5), Zen(:,5), dist(:,5), time(:,5)];

for i = 1:5
    min_dist(1,i) = min(dist(:,i));
    min_dist_row(1,i) = find(dist(:,i) == min_dist(1,i));
    min_dist_time(1,i) = time (min_dist_row(1,i),i);
    min_dist_zen(1,i) = Zen (min_dist_row(1,i),i);
end

for i = 1:5
    min_zen(1,i) = min(Zen(:,i));
    min_zen_row(1,i) = find(Zen(:,i) == min_zen(1,i));
    min_zen_time(1,i) = time (min_zen_row(1,i),i);
end

min_zen=min_zen*180/pi;
min_zen_time=min_zen_time/3600;

min_dist_time=min_dist_time/3600;
min_dist_zen=min_dist_zen*180/pi;

figure;
plot(linspace(0,24,7200)' , azim(:,1))
grid on
title('Azimuth in time - ErSa II')
xlabel('Time (hour)')
ylabel('Azimuth')

figure;
plot(linspace(0,24,7200)' , azim(:,2))
grid on
title('Azimuth in time - GPS')
xlabel('Time (hour)')
ylabel('Azimuth')

figure;
plot(linspace(0,24,7200)' , azim(:,3))
grid on
title('Azimuth in time - CHAMP')
xlabel('Time (hour)')
ylabel('Azimuth')

figure;
plot(linspace(0,24,7200)' , azim(:,4))
grid on
title('Azimuth in time - Molniya')
xlabel('Time (hour)')
ylabel('Azimuth')

figure;
plot(linspace(0,24,7200)' , azim(:,5))
grid on
title('Azimuth in time - Beidou')
xlabel('Time (hour)')
ylabel('Azimuth')

 %% 
figure;
plot(linspace(0,24,7200)' , elev(:,1));
grid on
title('Elevation angle in time - Ersa II')
xlabel('Time (hour)')
ylabel('Elevation Angle')

figure;
plot(linspace(0,24,7200)' , elev(:,2));
grid on
title('Elevation angle in time - GPS')
xlabel('Time (hour)')
ylabel('Elevation Angle')

figure;
plot(linspace(0,24,7200)' , elev(:,3));
grid on
title('Elevation angle in time - CHAMP')
xlabel('Time (hour)')
ylabel('Elevation Angle')

figure;
plot(linspace(0,24,7200)' , elev(:,4));
grid on
title('Elevation angle in time - Molniya')
xlabel('Time (hour)')
ylabel('Elevation Angle')

figure;
plot(linspace(0,24,7200)' , elev(:,5));
grid on
title('Elevation angle in time - Beidou')
xlabel('Time (hour)')
ylabel('Elevation Angle')