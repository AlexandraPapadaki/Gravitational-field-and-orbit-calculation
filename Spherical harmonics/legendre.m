clc
clear all
close all
pnm=zeros(11,11);
angle=zeros(181,1);
pnm4_0=zeros(181,1);
pnm5_5=zeros(181,1);
pnm6_2=zeros(181,1);
pnm7_6=zeros(181,1);
 for i=2:181
    angle(i,1)=i-1;
 end
for q=1:181
 t=cos((q-1)*pi/180);
 a1=sqrt(3);
 %ypologismos kokkinis grammis
 pnm(1,1)=1;
 pnm(2,2)=a1*sin((q-1)*pi/180)*pnm(1,1);
 for k=2:10
    a=sqrt((2*k+1)/(2*k));
    pnm(k+1,k+1)=a*sin((q-1)*pi/180)*pnm(k,k);
 end
 %ypologismos prasinis grammis
 for k=2:11
    n=k-1;
    m=k-2;
    b=sqrt(((2*n+1)*(2*n-1))/((n+m)*(n-m)));
    pnm(k,k-1)=b*t*pnm(k-1,k-1);
 end
 %ypologismos mple grammis
 for k=1:9
    for j=k+2:11
        n=j-1;
        m=k-1;
        b=sqrt(((2*n+1)*(2*n-1))/((n+m)*(n-m)));
        c=sqrt(((2*n+1)*(n-m-1)*(n+m-1))/((2*n-3)*(n+m)*(n-m)));
        pnm(j,k)=b*t*pnm(j-1,k)-c*pnm(j-2,k);    
    end
 end
 pnm4_0(q,1)=pnm(5,1);
 pnm5_5(q,1)=pnm(6,6);
 pnm6_2(q,1)=pnm(7,3);
 pnm7_6(q,1)=pnm(8,7);
end
figure
plot(angle,pnm4_0,'r');
% xlabel('colatitude è (degrees)')
% ylabel('m=0 and n=4')
grid on
% figure
hold on
plot(angle,pnm5_5,'b');
% xlabel('colatitude è (degrees)')
% ylabel('m=5 and n=5')
grid on
% figure
hold on
plot(angle,pnm6_2,'g');
% xlabel('colatitude è (degrees)')
% ylabel('m=2 and n=6')
grid on
% figure
hold on
plot(angle,pnm7_6,'y');
xlabel('colatitude è (degrees)')
ylabel('Pnm')
grid on
legend('m=0 and n=4','m=5 and n=5','m=2 and n=6','m=6 and n=7')