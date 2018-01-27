clc
clear all
close all

angle=zeros(181,1);
Pnm4_0=zeros(181,1);
Pnm5_5=zeros(181,1);
Pnm6_2=zeros(181,1);
Pnm7_6=zeros(181,1);

for i=2:181
    angle(i,1)=i-1;
end
 
for theta=0:180
    n = 4;
    m = 0;
    Pnm4_0(theta+1,1) = legendre(n,m,theta);
    
    n = 5;
    m = 5;
    Pnm5_5(theta+1,1) = legendre(n,m,theta);
    
    n = 6;
    m = 2;
    Pnm6_2(theta+1,1) = legendre(n,m,theta);
    
    n = 7;
    m = 6;
    Pnm7_6(theta+1,1) = legendre(n,m,theta);
end

figure
plot(angle,Pnm4_0,'r');
% xlabel('colatitude è (degrees)')
% ylabel('m=0 and n=4')
grid on
% figure
hold on
plot(angle,Pnm5_5,'b');
% xlabel('colatitude è (degrees)')
% ylabel('m=5 and n=5')
grid on
% figure
hold on
plot(angle,Pnm6_2,'g');
% xlabel('colatitude è (degrees)')
% ylabel('m=2 and n=6')
grid on
% figure
hold on
plot(angle,Pnm7_6,'y');
xlabel('colatitude è (degrees)')
ylabel('Pnm')
grid on
legend('m=0 and n=4','m=5 and n=5','m=2 and n=6','m=6 and n=7')