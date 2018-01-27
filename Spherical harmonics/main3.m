clear all
clc

r = 6378;
R = 6378;

GM = 398600.44;

T0 = GM/R;
T1 = 0;

%% n = 2
DC2(1,1) = 0.002*10^(-6);
DC2(1,2) = 0;
DC2(1,3) = 2.439*10^(-6);

DS2(1,1) = 1;
DS2(1,2) = 0;
DS2(1,3) = -1.4*10^(-6);
n2 = 2;
T2 = zeros(91,91);
for theta = 0:2:180
    for lamda = 0:4:360
        a = 0;
        for m2 = 0:n2
            Pnm = legendre(n2, m2, theta);
            a = a + Pnm*(DC2(1,m2+1)*cos(m2*lamda) + DS2(1,m2+1)*sin(m2*lamda));
        end
        T2(theta/2+1 , lamda/4+1) = GM/R*1*a;
    end
end

%% n = 3
DC3(1,1) = 0.957*10^(-6);
DC3(1,2) = 2.029*10^(-6);
DC3(1,3) = 0.904*10^(-6);
DC3(1,4) = 0.721*10^(-6);

DS3(1,1) = 1;
DS3(1,2) = 0.249*10^(-6);
DS3(1,3) = -0.619*10^(-6);
DS3(1,4) = 1.414*10^(-6);

n3 = 3;
T3 = zeros(91,91);
for theta = 0:2:180
    for lamda = 0:4:360
        a1 = 0;
        for m3 = 0:n3
            Pnm = legendre(n3, m3, theta);
            a1 = a1 + Pnm*(DC3(1,m3+1)*cos(m3*lamda) + DS3(1,m3+1)*sin(m3*lamda));
        end
        T3(theta/2+1 , lamda/4+1) = GM/R*1*a1;
    end
end


