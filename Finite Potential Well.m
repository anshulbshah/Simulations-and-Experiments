%% 
clear all
close all
V0 = 2*1.6e-19;
L = 5*1e-9;
radius = sqrt(2*9.1e-31*V0*L^2*(2*pi)^2/6.626e-34^2);
circle(0,0,radius);
kL = 0:0.00005:36;
KapL_even = kL.*tan(kL);
ev = find(abs(KapL_even-sqrt((radius^2)-kL.^2)) <= 1e-2);
ev2 = find(abs(KapL_even) <= radius);
KapL_even2 = KapL_even(ev);  
KapL_even = KapL_even(ev2);  
KapL_odd = -1.*kL.*cot(kL);
od = find(abs(KapL_odd-sqrt((radius^2)-kL.^2)) <= 1e-2);
od2 = find(abs(KapL_odd) <= radius);
KapL_odd2 = KapL_odd(od); 
KapL_odd = KapL_odd(od2); 
figure(1)
axis([0 40 0 40])
plot(kL(ev),KapL_even2,'ro');
hold on;
plot(kL(od),KapL_odd2,'bo');
hold on;
plot(kL(ev2),KapL_even,'--');
hold on;
plot(kL(od2),KapL_odd,'--');
xlabel('k L','FontSize',15);
ylabel('\kappa L','FontSize',15);
%% 
%State 1
%Even Parity
solutions_even = [kL(ev(1:5));KapL_even2(1:5)];
A_by_C = exp(-solutions_even(2,:))./cos(solutions_even(1,:));
C = (1/sqrt(L))*((A_by_C.^2).*(sin(2.*solutions_even(1,:))./(2*solutions_even(1,:))+1) + exp(-2.*solutions_even(2,:))/(solutions_even(2,:))).^(-2)

solutions_odd = [kL(od(1:5));KapL_odd2(1:5)];
A_by_C2 = exp(-solutions_odd(2,:))./sin(solutions_odd(1,:));
C2 = (1/L)*((A_by_C2.^2).*(-L/4*sin(2.*solutions_odd(1,:))./solutions_odd(1,:)+L) + exp(-2.*solutions_odd(2,:))/(2.*solutions_odd(2,:))).^(-2)

x_L = -1.5*L:0.01*L:-1*L;
x_R = L:0.01*L:1.5*L;
x_M = -L:0.01*L:L;
plot(x_M,A_by_C(1)*C(1)*cos(solutions_even(1,1).*x_M/L));
hold on;
plot(x_L,C(1)*exp(solutions_even(2,1).*x_L/L));
hold on;
plot(x_R,C(1)*exp(-solutions_even(2,1).*x_R/L));

figure;
plot(x_M,A_by_C2(1)*C2(1)*sin(solutions_odd(1,1).*x_M/L)*L);
hold on;
plot(x_L,-C2(1)*exp(solutions_odd(2,1).*x_L/L)*L);
hold on;
plot(x_R,C2(1)*exp(-solutions_odd(2,1).*x_R/L)*L);

figure;
plot(x_M,-1*A_by_C(2)*C(2)*1e-10*cos(solutions_even(1,2).*x_M/L)*L);
hold on;
plot(x_L,-C(2)*1e-10*exp(solutions_even(2,2).*x_L/L)*L);
hold on;
plot(x_R,-C(2)*1e-10*exp(-solutions_even(2,2).*x_R/L)*L);