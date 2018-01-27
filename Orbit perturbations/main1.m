clear all
close all
clc

%% 2i
iter_a1 = (31412500-6838000)/1500;
dw1 = zeros(iter_a1,1);
dW1 = zeros(iter_a1,1);
dM1 = zeros(iter_a1,1);
a = zeros(iter_a1,1);
n1 = zeros(iter_a1,1);

a1 = 6838000;
for k = 1:iter_a1
    a(k,1) = a1;
    e1 = 0.01;
    i1 = 55*pi()/180;
    [dw1(k,1), dW1(k,1), dM1(k,1), n1(k,1)] = disturb(a1,e1,i1);
    a1 = a1 + 1500;
end
plot(a,dw1,'r')
hold on
ylabel('dw/dt');
xlabel('a');
grid on
title('Quantity dw/dt as a function of a');
figure

plot(a,dW1,'r')
hold on
ylabel('dW/dt');
xlabel('a');
grid on
title('Quantity dW/dt as a function of a');
figure

plot(a, dM1-n1,'r')
hold on
ylabel('dM/dt-n');
xlabel('a');
grid on
title('Quantity dM/dt-n as a function of a');

%% 2ii
iter_e2 = (0.9-0.00001)/0.0001;
dw2 = zeros(floor(iter_e2),1);
dW2 = zeros(floor(iter_e2),1);
dM2 = zeros(floor(iter_e2),1);
e = zeros(iter_e2,1);
n2 = zeros(iter_e2,1);

e2 = 0.00001;
for k = 1:floor(iter_e2)
    e(k,1) = e2;
    a2 = 6838000;
    i2 = 55*pi()/180;
    [dw2(k,1), dW2(k,1), dM2(k,1), n2(k,1)] = disturb(a2,e2,i2);
    e2 = e2 + 0.0001;
end

figure
plot(e,dw2,'r')
hold on
ylabel('dw/dt');
xlabel('e');
grid on
title('Quantity dw/dt as a function of e');
figure

plot(e,dW2,'r')
hold on
ylabel('dW/dt');
xlabel('e');
grid on
title('Quantity dW/dt as a function of e');
figure

plot(e, dM2-n2,'r')
hold on
ylabel('dM/dt-n');
xlabel('e');
grid on
title('Quantity dM/dt-n as a function of e');

%% 2iii
iter_i3 = (90-0)/0.01;
dw3 = zeros(iter_i3,1);
dW3 = zeros(iter_i3,1);
dM3 = zeros(iter_i3,1);
i = zeros(iter_i3,1);
n3 = zeros(iter_i3,1);

i3 = 0;
for k = 1:iter_i3
    i(k,1) = i3*pi()/180;
    a3 = 6838000;
    e3 =  0.01;
    [dw3(k,1), dW3(k,1), dM3(k,1), n3(k,1)] = disturb(a3,e3,i3);
    i3 = i3 + 0.01*pi()/180;
end

figure
plot(i*180/pi(),dw3,'r')
hold on
ylabel('dw/dt');
xlabel('i');
grid on
title('Quantity dw/dt as a function of i');
figure

plot(i*180/pi(),dW3,'r')
hold on
ylabel('dW/dt');
xlabel('i');
grid on
title('Quantity dW/dt as a function of i');
figure

plot(i*180/pi(), dM3-n3,'r')
hold on
ylabel('dM/dt-n');
xlabel('i');
grid on
title('Quantity dM/dt-n as a function of i');