clear;clc
options=odeset('RelTol',1e-13,'AbsTol',1e-14,'Refine',8); 
%Number of points to evaluate
N = 10000;
%integration interval
tspan = linspace(0,5,N); %when ga = 1, mu = 2
% tspan = linspace(0,4*20,N); %when ga = 0.1, mu = 4

% initial states
spin_up = [1;0];
spin_down = [0;1];

ga = 0.1;
mu = 4;
%solve the first order schrodinger equation
% [time1, state1] = ode45(@(t,psi) H1(t,psi,ga,mu),tspan,spin_up);
% [time2, state2] = ode45(@(t,psi) H1(t,psi,ga,mu),tspan,spin_down);

[time1, state1] = ode45(@(t,psi) Hermitian_H(t,psi,ga,mu,5),tspan,spin_up);
[time2, state2] = ode45(@(t,psi) Hermitian_H(t,psi,ga,mu,5),tspan,spin_down);



%array of population of spin up/down 
spinuparray = abs(state1(:,1)).^2; %project spinup and mod square
spindownarray = abs(state1(:,2)).^2;

plot(time1,spinuparray,'b-',time2,spindownarray,'r--','linewidth',2)
grid on
set(gca,'FontSize',50)
xlabel('t','FontSize',50)
ylabel('Population','FontSize',50)
% axis([0 5 -0.1 1.7]) %when ga=1, mu=2
% axis([0 20 0 4]) %when ga=0.1, mu=4