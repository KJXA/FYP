clear;clc
%Number of points to evaluate
N = 1000;
%integration interval
tspan = linspace(0,1,N); 
% initial states
spin_up = [1;0];
spin_down = [0;1];

ga = 0.1;
mu = 4;
%solve the first order schrodinger equation
[time1, state1] = ode45(@(t,psi) H1(t,psi,ga,mu),tspan,spin_up);
[time2, state2] = ode45(@(t,psi) H1(t,psi,ga,mu),tspan,spin_down);

%solve quadratic eigenvalue equation
y11 = state1(:,1);y21 = state1(:,2);
y12 = state2(:,1);y22 = state2(:,2);
b = -y11-y22;
c = y11.*y22-y12.*y21;
eigen1 = (-b+sqrt(b.^2-4.*c))/2;
eigen2 = (-b-sqrt(b.^2-4.*c))/2;
real_eigen1 = real(eigen1);
real_eigen2 = real(eigen2);
plot(time1,real_eigen1,'b-',time2,real_eigen2,'r--','linewidth',2)
grid on
set(gca,'FontSize',50)
xlabel('\textit{t}', 'Interpreter','Latex','FontSize',50)
ylabel('$Re(\lambda_{U})$ ', 'Interpreter','Latex','FontSize',50)
% axis([0 1 0 2]) % for ga=1 mu=2
% axis([0 5 0 4]) % for ga=0.1, mu=4