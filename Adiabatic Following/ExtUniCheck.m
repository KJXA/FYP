clear;clc;
%Number of points to evaluate
N = 1000;
T=100;
ga=5;
mu=1.2;
%integration interval
tspan = linspace(0,T,N); 

% initial states
spin_up = [1;0];
spin_down = [0;1];
options=odeset('RelTol',1e-13,'AbsTol',1e-14,'Refine',8);
%get cyclic states
% [time1, state1]=ode45(@(t,y) H2(t,y,mu,T),tspan,spin_up);
% [time2, state2]=ode45(@(t,y) H2(t,y,mu,T),tspan,spin_down);
[time1, state1]=ode45(@(t,y) H2(t,y,ga,mu,T),tspan,spin_up,options);
[time2, state2]=ode45(@(t,y) H2(t,y,ga,mu,T),tspan,spin_down,options);

u11=state1(end,1);
u12=state2(end,1);
u21=state1(end,2);
u22=state2(end,2);

[vec,values] = eig([u11,u12;u21,u22], 'vector');
abs(values).^2

%Extended Unitarity still exists for T=100, mu > 0.33 i.e hopping occurs
%But what about converse? No hopping even when no extended unitarity.

%get instantaneous eigenvalues
inenergy1 = zeros(1,N);
inenergy2 = zeros(1,N);
index = 1;
for t = tspan
[invec,invalue] = eig(h2_matrix(t,mu,T),'vector');
inenergy1(:,index)=invalue(1);
inenergy2(:,index)=invalue(2);
index = index + 1;
end
ax1=subplot(1,2,1);
plot(ax1,tspan,real(inenergy1),'r',tspan,real(inenergy2),'b')
xlabel('t');
ylabel('Eigenenergy');
ax2=subplot(1,2,2);
plot(ax2,tspan,imag(inenergy1),'r',tspan,imag(inenergy2),'b')
xlabel('t');
ylabel('Eigenenergy');