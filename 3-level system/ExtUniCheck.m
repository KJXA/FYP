clear;clc
N = 1000;
T=100;
ga=15;
mu=0.8;

%integration interval
tspan = linspace(0,T,N); 

% initial states
spin_up = [1;0;0];
spin_mid = [0;1;0];
spin_down = [0;0;1];
options=odeset('RelTol',1e-15,'AbsTol',1e-20,'Refine',8);
% get cyclic states
[time1, state1] = ode45(@(t,psi) H_1(t,psi,mu,T),tspan,spin_up,options);
[time2, state2] = ode45(@(t,psi) H_1(t,psi,mu,T),tspan,spin_mid,options);
[time3, state3] = ode45(@(t,psi) H_1(t,psi,mu,T),tspan,spin_down,options);

%get eigenvalues
temp1=transpose(state1(end,:));
temp2=transpose(state2(end,:));
temp3=transpose(state3(end,:));
floquet = [temp1 temp2 temp3];
[vec,eigen] = eig(floquet,'vector');

abs(eigen).^2

%get instantaneous eigenvalues
inenergy1 = zeros(1,N);
inenergy2 = zeros(1,N);
inenergy3 = zeros(1,N);
index = 1;
for t = tspan
[invec,invalue] = eig(h1(t,mu,T),'vector');
inenergy1(:,index)=invalue(1);
inenergy2(:,index)=invalue(2);
inenergy3(:,index)=invalue(3);
index = index + 1;
end
ax1=subplot(1,2,1);
plot(ax1,tspan,real(inenergy1),'r',tspan,real(inenergy2),'b',tspan,real(inenergy3),'g')
xlabel('t');
ylabel('Eigenenergy');
ax2=subplot(1,2,2);
plot(ax2,tspan,imag(inenergy1),'r',tspan,imag(inenergy2),'b',tspan,imag(inenergy3),'g')
xlabel('t');
ylabel('Eigenenergy');