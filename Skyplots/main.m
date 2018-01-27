close all
clear all
clc

%% right ascension hours--> degrees --> rad
a1 = 05 + 8/60 + 42.36345199/3600;
a1 = a1*15*pi()/180;

a2 = 11 + 13/60 + 58.69508359/3600;
a2 = a2*15*pi()/180;

a3 = 17 + 39/60 + 27.39049252/3600;
a3 = a3*15*pi()/180;

a4 = 11 + 03/60 + 52.22167171/3600;
a4 = a4*15*pi()/180;

%% declination in degrees --> rad
d1 = (84 + 32/60 + 04.5440155/3600)*pi()/180;

d2 = (14 + 42/60 + 26.9525965/3600)*pi()/180;

d3 = (49 + 55/60 + 03.3684410/3600)*pi()/180;

d4 = (-53 + 57/60 + 00.6966293/3600)*pi()/180;

coord = [a1, d1; a2, d2; a3, d3; a4, d4];

%% calculate x(gcrs)
x = zeros(4,3);
for i = 1 :4
    x(i,1) = cos(coord(i,2))*cos(coord(i,1));
    x(i,2) = cos(coord(i,2))*sin(coord(i,1));
    x(i,3) = sin(coord(i,2));
end

phi = 51+30/60+01/3600;
phi = phi*pi()/180;

lamda = 13+24/60+29/3600;
lamda = lamda*pi()/180;

coord = x;

Y =2015;
M = 11;
D = 2;

M1 = [1 0 0; 0 -1 0; 0 0 1];

for i = 1:4
    %% calculate N, P
    dt = 15/60;
    DUT1 = 0.1/3600;
    CET = 0;
    UT1 = CET - 1 + DUT1;
    for j = 1:97
        a = Y + floor((M + 9)/12);
        JD = 367*Y - floor(7*a/4) + floor(275*M/9)+ D + 1721014 + UT1/24 - 0.5;
        d = JD - 2451545.0;
    
        T = d/36525;
        
        e =(84381.448/3600-46.8150/3600*T)*pi/180;

        f1 = (125 - 0.05295*d)*pi()/180;
        f2 = (200.9 + 1.97129*d)*pi()/180;
        dy = (-0.0048*sin(f1) - 0.0004*sin(f2))*pi()/180;
        de = (0.0026*cos(f1) - 0.0002*cos(f2))*pi()/180;

        N = R1(-e)*R3(dy)*R1(de+e);

        tha = ((2004.3109*T - 0.42665*T^2)/3600)*pi()/180;
        za = ((2306.2181*T + 1.09468*T^2)/3600)*pi()/180;
        zita = ((2306.2181*T + 0.30188*T^2)/3600)*pi()/180;
        P = R3(-za)*R2(tha)*R3(-zita);

        %% Time scales
        if j < 90
            n = -16;
            GMST = (UT1*3600 + 24110.54841 + 8640184.812866*T + 0.093104*T^2 - 6.2*10^(-6)*T^3)/3600 + n*24;
        else 
            n = -17;
            GMST = (UT1*3600 + 24110.54841 + 8640184.812866*T + 0.093104*T^2 - 6.2*10^(-6)*T^3)/3600 + n*24;
        end
       
        GMST_out(j,1) = GMST;
        
        GAST = GMST*pi()/12 + (dy*cos(e+de))/15;
        %%
        x = M1*R2(90*pi()/180-phi)*R3(lamda)*R3(GAST)*N*P*coord(i,:)';
        xlocal (j,3*(i-1)+1) = x (1,1);
        xlocal (j,3*(i-1)+2) = x (2,1);
        xlocal (j,3*(i-1)+3) = x (3,1);
        
        UT1 = UT1 + dt;
    end
end

%% calculation of elevation e, azimuth A
%% source 1
for i = 1:97
    c = acot(abs(xlocal(i,1))/abs(xlocal(i,2)));
    if xlocal(i,1)>0 && xlocal(i,2)>0
        eA1(i,2) = c;
    else if xlocal(i,1)<0 && xlocal(i,2)>0
            eA1(i,2) = pi()-c;
        else if xlocal(i,1)<0 && xlocal(i,2)<0
                eA1(i,2) = pi()+c;
            else if xlocal(i,1)>0 && xlocal(i,2)<0
                   eA1(i,2) = 2*pi()-c;
                end
            end
        end
    end
    eA1(i,1) = asin(xlocal(i,3));
end

%% source 2
for i = 1:97
    c = acot(abs(xlocal(i,4))/abs(xlocal(i,5)));
    if xlocal(i,4)>0 && xlocal(i,5)>0
        eA2(i,2) = c;
    else if xlocal(i,4)<0 && xlocal(i,5)>0
            eA2(i,2) = pi()-c;
        else if xlocal(i,4)<0 && xlocal(i,5)<0
                eA2(i,2) = pi()+c;
            else if xlocal(i,4)>0 && xlocal(i,5)<0
                   eA2(i,2) = 2*pi()-c;
                end
            end
        end
    end
    eA2(i,1) = asin(xlocal(i,6));
end

%% source 3
for i = 1:97
    c = acot(abs(xlocal(i,7))/abs(xlocal(i,8)));
    if xlocal(i,7)>0 && xlocal(i,8)>0
        eA3(i,2) = c;
    else if xlocal(i,7)<0 && xlocal(i,8)>0
            eA3(i,2) = pi()-c;
        else if xlocal(i,7)<0 && xlocal(i,8)<0
                eA3(i,2) = pi()+c;
            else if xlocal(i,7)>0 && xlocal(i,8)<0
                   eA3(i,2) = 2*pi()-c;
                end
            end
        end
    end
    eA3(i,1) = asin(xlocal(i,9));
end

%% source 4
for i = 1:97
    c = acot(abs(xlocal(i,10))/abs(xlocal(i,11)));
    if xlocal(i,10)>0 && xlocal(i,11)>0
        eA4(i,2) = c;
    else if xlocal(i,10)<0 && xlocal(i,11)>0
            eA4(i,2) = pi()-c;
        else if xlocal(i,10)<0 && xlocal(i,11)<0
                eA4(i,2) = pi()+c;
            else if xlocal(i,10)>0 && xlocal(i,11)<0
                   eA4(i,2) = 2*pi()-c;
                end
            end
        end
    end
    eA4(i,1) = asin(xlocal(i,12));
end

%Diurnal paths
plot(1:97,eA1(:,1),'r');
hold on;
plot(1:97,eA2(:,1),'b');
hold on;
plot(1:97,eA3(:,1),'g');
hold on;
plot(1:97,eA4(:,1),'y');
grid on
legend('ICRF J050842.3 + 843204', 'ICRF J111358.6 + 144226', 'ICRF J173927.3 + 495503', 'ICRF J110352.2 + 535700')

%other plot
figure
x1 = xlocal(:,1);
y1 = xlocal(:,2);
z1 = xlocal(:,3);
plot3(x1,y1,z1,'.r');
hold on
x2 = xlocal(:,4);
y2 = xlocal(:,5);
z2 = xlocal(:,6);
plot3(x2,y2,z2,'.b');
hold on
x3 = xlocal(:,7);
y3 = xlocal(:,8);
z3 = xlocal(:,9);
plot3(x3,y3,z3,'.g');
hold on
x4 = xlocal(:,10);
y4 = xlocal(:,11);
z4 = xlocal(:,12);
plot3(x4,y4,z4,'.y');
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
hold on
grid on 
legend('ICRF J050842.3 + 843204', 'ICRF J111358.6 + 144226', 'ICRF J173927.3 + 495503', 'ICRF J110352.2 + 535700')

%% skyplots
figure;
skyplot(eA1(:,2)*180/pi(),eA1(:,1)*180/pi(),'.r')
hold on

skyplot(eA2(:,2)*180/pi(),eA2(:,1)*180/pi(),'.b')
hold on

skyplot(eA3(:,2)*180/pi(),eA3(:,1)*180/pi(),'.g')
hold on

skyplot(eA4(:,2)*180/pi(),eA4(:,1)*180/pi(),'.y')
hold on 

legend('ICRF J050842.3 + 843204', 'ICRF J111358.6 + 144226', 'ICRF J173927.3 + 495503', 'ICRF J110352.2 + 535700', 'Celestial pole')
