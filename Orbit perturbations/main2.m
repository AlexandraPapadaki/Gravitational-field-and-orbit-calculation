clear all
close all
clc

a_champ = 6838000;
e_champ = 0.004;
i_champ = 87*pi()/180;
w_champ = 0;
W_champ = 0;
M_champ = 0;

GM = 398600.44;
w_earth=2*pi()/86164;

initial = [a_champ e_champ i_champ w_champ W_champ M_champ];

%% Perturbed Keplerian elements
[dw, dW, dM, n] = disturb(a_champ,e_champ,i_champ);

T = 4*pi()/(2*n);

dt = 2*T/720;

a_per = zeros(721,6);
a_per(:,1) = a_champ;
a_per(:,2) = e_champ;
a_per(:,3) = i_champ;
Dt = zeros(721,6);
for k = 1:721
    M_champ=((k-1)*pi/180)';
    Dt(k,1) = (k-1)*dt;
    
    a_per(k,4) = w_champ + dw*Dt(k,1);
    a_per(k,5) = W_champ + dW*Dt(k,1);
    a_per(k,6) = M_champ + dM*Dt(k,1);
end

%% Perturbed orbit in Space-fixed
E_f=zeros(721,1);
r_pos=zeros(721,3);
Eo(:,1) = a_per(:,6);
E(:,1) = a_per(:,6) + e_champ*sin(Eo(:,1));
    for i = 1:50
        E(:,i+1) = a_per(:,6) + e_champ*sin(E(:,i));
        if abs(E(:,i+1)-E(:,i))<10^(-8)
        break
        end
    end
E_f(:,1) = E(:,end);

for k = 1:721
    K(:,1)=[cos(a_per(k,5));sin(a_per(k,5));0];
    C0(:,1)=[sin(a_per(k,5))*sin(a_per(k,3));-cos(a_per(k,5))*sin(a_per(k,3));cos(a_per(k,3))];
    P(:,1)=cos(a_per(k,4))*K(:,1)+sin(a_per(k,4))*(cross(C0(:,1),K(:,1)));
    Q(:,1)=-sin(a_per(k,4))*K(:,1)+cos(a_per(k,4))*(cross(C0(:,1),K(:,1)));
    
    f_rad=2*atan(tan(E_f(k,1)/2)*sqrt((1+a_per(k,2))/(1-a_per(k,2))));
    r=a_per(k,1)*(1-a_per(k,2)*cos(E_f(k,1)));
    p=r*(1+a_per(k,2)*cos(f_rad));
    
    r_pos(k,1:3)=r*(cos(f_rad)*P(:,1)+sin(f_rad)*Q(:,1));
end

%% Unperturbed orbit in Space-fixed
a_unper=zeros(721,6);
for k=1:721
    a_unper(k,1:5) = initial(1:5);
    a_unper(k,6)=(k-1)*pi()/180;
end

E_f_unper=zeros(721,1);
r_pos_un=zeros(721,3);
Eo_un(:,1) = a_unper(:,6);
E_un(:,1) = a_unper(:,6) + e_champ*sin(Eo_un(:,1));
    for i = 1:50
        E_un(:,i+1) = a_unper(:,6) + e_champ*sin(E_un(:,i));
        if abs(E_un(:,i+1)-E_un(:,i))<10^(-8)
        break
        end
    end
E_f_unper(:,1) = E_un(:,end);

for k = 1:721
    K(:,1)=[cos(a_unper(k,5));sin(a_unper(k,5));0];
    C0(:,1)=[sin(a_unper(k,5))*sin(a_unper(k,3));-cos(a_unper(k,5))*sin(a_unper(k,3));cos(a_unper(k,3))];
    P(:,1)=cos(a_unper(k,4))*K(:,1)+sin(a_unper(k,4))*(cross(C0(:,1),K(:,1)));
    Q(:,1)=-sin(a_unper(k,4))*K(:,1)+cos(a_unper(k,4))*(cross(C0(:,1),K(:,1)));
    
    f_rad_unper=2*atan(tan(E_f_unper(k,1)/2)*sqrt((1+a_unper(k,2))/(1-a_unper(k,2))));
    r=a_unper(k,1)*(1-a_unper(k,2)*cos(E_f_unper(k,1)));
    p=r*(1+a_unper(k,2)*cos(f_rad_unper));
    
    r_pos_un(k,1:3)=r*(cos(f_rad_unper)*P(:,1)+sin(f_rad_unper)*Q(:,1));
end

%% Perturbed and Unperturbed orbit in Earth-fixed
E_int=1;
pos_true = zeros(721,3);
pos_true_un = zeros(721,3);
for k = 1:721
    t_dE=(k-1)*T*E_int/360;
    dw_ang=w_earth*t_dE;
    R_3=[cos(-dw_ang) -sin(-dw_ang) 0;sin(-dw_ang) cos(-dw_ang) 0;0 0 1];
    pos_true(k,1:3)=(r_pos(k,1:3)*R_3');
    pos_true_un(k,1:3)=(r_pos_un(k,1:3)*R_3');
end

lat=zeros(721,1);
lon=zeros(721,1);
h=zeros(721,1);
lat_un=zeros(721,1);
lon_un=zeros(721,1);
h_un=zeros(721,1);
for i = 1:721
    [lat(i,1),lon(i,1),h(i,1)] = xyz2blh(pos_true(i,1),pos_true(i,2),pos_true(i,3));
    [lat_un(i,1),lon_un(i,1),h_un(i,1)] = xyz2blh(pos_true_un(i,1),pos_true_un(i,2),pos_true_un(i,3));
end

%% Perturbed and Unperturbed orbit in Topocentric
X_ber=3783.26649;
Y_ber=901.6469;
Z_ber=5038.24814;

l_local=atan2(Y_ber,X_ber);
f_local=atan(Z_ber/sqrt(X_ber^2+Y_ber^2));
rotations=[-1,0,0;0,1,0;0,0,1]*[cos(pi/2-f_local),0,sin(pi/2-f_local);0,1,0;-sin(pi/2-f_local),0,cos(pi/2-f_local)]*[cos(l_local),-sin(l_local),0;sin(l_local),cos(l_local),0;0,0,1];

Dr = zeros(721,1);
Dr_un = zeros(721,1);
for k=1:721
Dr(k,1)=(pos_true(k,1)/10^3)-X_ber;
Dr(k,2)=(pos_true(k,2)/10^3)-Y_ber;
Dr(k,3)=(pos_true(k,3)/10^3)-Z_ber;
Dr_un(k,1)=(pos_true_un(k,1)/10^3)-X_ber;
Dr_un(k,2)=(pos_true_un(k,2)/10^3)-Y_ber;
Dr_un(k,3)=(pos_true_un(k,3)/10^3)-Z_ber;
end

El = zeros(721,1);
El_un = zeros(721,1);
Az = zeros(721,1);
Az_un = zeros(721,1);
pos_cel=zeros(721,3);
pos_cel_un=zeros(721,3);
v=zeros(720,3);
v_un=zeros(720,3);
for k=1:721
pos_cel(k,:)=(rotations*Dr(k,:)')';
Az(k,1)=atan2(pos_cel(k,2),pos_cel(k,1));
El(k,1)=atan2(pos_cel(k,3),sqrt(pos_cel(k,1)^2+pos_cel(k,2)^2));

pos_cel_un(k,:)=(rotations*Dr_un(k,:)')';
Az_un(k,1)=atan2(pos_cel_un(k,2),pos_cel_un(k,1));
El_un(k,1)=atan2(pos_cel_un(k,3),sqrt(pos_cel_un(k,1)^2+pos_cel_un(k,2)^2));
end
pos_cel_true=[pos_cel(:,1)+X_ber,pos_cel(:,2)+Y_ber,pos_cel(:,3)+Z_ber];
pos_cel_true_un=[pos_cel_un(:,1)+X_ber,pos_cel_un(:,2)+Y_ber,pos_cel_un(:,3)+Z_ber];

for k=1:720
v(k,1)=(pos_cel(k+1,1)-pos_cel(k,1))/dt;
v(k,2)=(pos_cel(k+1,2)-pos_cel(k,2))/dt;
v(k,3)=(pos_cel(k+1,3)-pos_cel(k,3))/dt;

v_un(k,1)=(pos_cel_un(k+1,1)-pos_cel_un(k,1))/dt;
v_un(k,2)=(pos_cel_un(k+1,2)-pos_cel_un(k,2))/dt;
v_un(k,3)=(pos_cel_un(k+1,3)-pos_cel_un(k,3))/dt;
end

Az_deg=Az*180/pi;
Az_deg_un=Az_un*180/pi;
El_deg=El*180/pi;
El_deg_un=El_un*180/pi;

%% Plots
plot(Dt,r_pos(:,1)-r_pos_un(:,1));
xlabel('time(sec)')
ylabel('X un(space fixed-m)-X (space fixed-m)')
grid on 
figure
plot(Dt,r_pos(:,2)-r_pos_un(:,2));
xlabel('time(sec)')
ylabel('Y un(space fixed-m)-Y (space fixed-m)')
grid on 
figure
plot(Dt,r_pos(:,3)-r_pos_un(:,3));
xlabel('time(sec)')
ylabel('Z un(space fixed-m)-Z (space fixed-m)')
grid on 
figure
plot(Dt,lon_un-lon);
xlabel('time(sec)')
ylabel('lon un(deg,m)-lon(deg,m)')
grid on 
figure
plot(Dt,lat_un-lat);
xlabel('time(sec)')
ylabel('lat un(deg,m)-lat(deg,m)')
grid on 
figure
plot(Dt,pos_cel_true_un(:,1)-pos_cel_true(:,1));
xlabel('time(sec)')
ylabel('rx un (km)-rx (km)')
grid on 
figure
plot(Dt,pos_cel_true_un(:,2)-pos_cel_true(:,2));
xlabel('time(sec)')
ylabel('ry un (km)-ry (km)')
grid on 
figure
plot(Dt,pos_cel_true_un(:,3)-pos_cel_true(:,3));
xlabel('time(sec)')
ylabel('rz un (km)-rz (km)')
grid on 
figure
plot(Dt(2:721,1),v_un(:,1)-v(:,1));
xlabel('time(sec)')
ylabel('vx un (km/s)-vx (km/s)')
grid on 
figure
plot(Dt(2:721,1),v_un(:,2)-v(:,2));
xlabel('time(sec)')
ylabel('vy un (km/s)-vy (km/s)')
grid on 
figure
plot(Dt(2:721,1),v_un(:,3)-v(:,3));
xlabel('time(sec)')
ylabel('vz un (km/s)-vz (km/s)')
grid on 
figure
plot(Dt,Az_deg_un-Az_deg);
xlabel('time(sec)')
ylabel('Azimuth un (degrees)-Azimuth (degrees)')
grid on 
figure
plot(Dt,El_deg_un-El_deg);
xlabel('time(sec)')
ylabel('Elevation un (degrees)-Elevation (degrees)')
grid on 