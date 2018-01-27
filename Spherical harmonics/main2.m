clc
clear all
close all

Cnm4_0 = zeros(91,91);
Snm4_0 = zeros(91,91);

n = 4;
m = 0;
for theta = 0:2:180
    for lambda = 0:4:360
        [Cnm4_0(theta/2+1,lambda/4+1),Snm4_0(theta/2+1,lambda/4+1)]=fun(n,m,theta,lambda);
    end
end

Cnm5_5 = zeros(91,91);
Snm5_5 = zeros(91,91);

n = 5;
m = 5;
for theta = 0:2:180
    for lambda = 0:4:360
        [Cnm5_5(theta/2+1,lambda/4+1),Snm5_5(theta/2+1,lambda/4+1)]=fun(n,m,theta,lambda);
    end
end

Cnm7_6 = zeros(91,91);
Snm7_6 = zeros(91,91);

n = 7;
m = 6;
for theta = 0:2:180
    for lambda = 0:4:360
        [Cnm7_6(theta/2+1,lambda/4+1),Snm7_6(theta/2+1,lambda/4+1)]=fun(n,m,theta,lambda);
    end
end
 
% surf(Cnm4_0, Snm4_0);

theta = 0:2:180;
lambda = 0:4:360;

imagesc(lambda',theta',Cnm4_0);
hold on

figure 
imagesc(lambda',theta',Snm4_0);
hold on

figure
imagesc(lambda',theta',Cnm5_5);
hold on

figure 
imagesc(lambda',theta',Snm5_5);
hold on

figure
imagesc(lambda',theta',Cnm7_6);
hold on

figure 
imagesc(lambda',theta',Snm7_6);
hold on

figure 
sphere 
hold on

surf(lambda',theta',Cnm4_0);