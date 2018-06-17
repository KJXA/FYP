clear;clc;clf
%Number of points to evaluate
N = 10000;
T=100;
ga=1;
mu=0.7; %cylic states change position (T=500, mu = 0.2,0.4) (T=150, mu=0.2,***0.8***,1.2)
%integration interval
tspan = linspace(0,T,N);

% initial states
spin_up = [1;0;0];
spin_mid = [0;1;0];
spin_down = [0;0;1];
options=odeset('RelTol',1e-15,'AbsTol',1e-20,'Refine',8);
[time1, state1] = ode45(@(t,psi) H_1(t,psi,mu,T),tspan,spin_up,options);
[time2, state2] = ode45(@(t,psi) H_1(t,psi,mu,T),tspan,spin_mid,options);
[time3, state3] = ode45(@(t,psi) H_1(t,psi,mu,T),tspan,spin_down,options);

%array of cyclic states
cyclic_states1 = zeros(3,N);
cyclic_states2 = zeros(3,N);
cyclic_states3 = zeros(3,N);

temp1=transpose(state1(end,:));
temp2=transpose(state2(end,:));
temp3=transpose(state3(end,:));
floquet = [temp1 temp2 temp3];
[init_cyc,eigval] = eig(floquet);

% initial cyclic states
init_cyc1=init_cyc(:,1); cyclic_states1(:,1)=init_cyc1;
init_cyc2=init_cyc(:,2); cyclic_states2(:,1)=init_cyc2;
init_cyc3=init_cyc(:,3); cyclic_states3(:,1)=init_cyc3;

for i = 2:N
    temp1=transpose(state1(i,:));
    temp2=transpose(state2(i,:));
    temp3=transpose(state3(i,:));
    floquet = [temp1 temp2 temp3];
    cyclic_states1(:,i)=floquet*init_cyc1;
    cyclic_states2(:,i)=floquet*init_cyc2;
    cyclic_states3(:,i)=floquet*init_cyc3;
end

%array of instantaneous eigenstates
invec1 = zeros(3,N);
invec2 = zeros(3,N);
invec3 = zeros(3,N);
index = 1;
for t = tspan
[invec,invalue] = eig(h1(t,mu,T));
invec1(:,index)=invec(:,1);
invec2(:,index)=invec(:,2);
invec3(:,index)=invec(:,3);
index = index + 1;
end
%get ratios
%%% b/a %%%
ratio_real_cy1 = real(cyclic_states1(2,:)./cyclic_states1(1,:)); ratio_real1 = real(invec1(2,:)./invec1(1,:));
ratio_im_cy1 = imag(cyclic_states1(2,:)./cyclic_states1(1,:)); ratio_im1 = imag(invec1(2,:)./invec1(1,:));
ratio_real_cy2 = real(cyclic_states2(2,:)./cyclic_states2(1,:)); ratio_real2 = real(invec2(2,:)./invec2(1,:));
ratio_im_cy2 = imag(cyclic_states2(2,:)./cyclic_states2(1,:)); ratio_im2 = imag(invec2(2,:)./invec2(1,:));
ratio_real_cy3 = real(cyclic_states3(2,:)./cyclic_states3(1,:)); ratio_real3 = real(invec3(2,:)./invec3(1,:));
ratio_im_cy3 = imag(cyclic_states3(2,:)./cyclic_states3(1,:)); ratio_im3 = imag(invec3(2,:)./invec3(1,:));

% %%% c/a %%%
% ratio_real_cy1 = real(cyclic_states1(2,:)./cyclic_states1(1,:)); ratio_real1 = real(invec1(2,:)./invec1(1,:));
% ratio_im_cy1 = imag(cyclic_states1(2,:)./cyclic_states1(1,:)); ratio_im1 = imag(invec1(2,:)./invec1(1,:));
% ratio_real_cy2 = real(cyclic_states2(2,:)./cyclic_states2(1,:)); ratio_real2 = real(invec2(2,:)./invec2(1,:));
% ratio_im_cy2 = imag(cyclic_states2(2,:)./cyclic_states2(1,:)); ratio_im2 = imag(invec2(2,:)./invec2(1,:));
% ratio_real_cy3 = real(cyclic_states3(2,:)./cyclic_states3(1,:)); ratio_real3 = real(invec3(2,:)./invec3(1,:));
% ratio_im_cy3 = imag(cyclic_states3(2,:)./cyclic_states3(1,:)); ratio_im3 = imag(invec3(2,:)./invec3(1,:));
% %%% c/a %%%

%plot ratios
ax1=subplot(2,3,1);
plot(ax1,tspan,ratio_real_cy1,'b-',tspan,ratio_real1,'r--',tspan,ratio_real2,'g--',tspan,ratio_real3,'m--') %Re(cyclic1) b/a
grid on;grid minor

ax2=subplot(2,3,2);
plot(ax2,tspan,ratio_real_cy2,'b-',tspan,ratio_real1,'r--',tspan,ratio_real2,'g--',tspan,ratio_real3,'m--') %Re(cyclic2) b/a
grid on;grid minor

ax3=subplot(2,3,3);
plot(ax3,tspan,ratio_real_cy3,'b-',tspan,ratio_real1,'r--',tspan,ratio_real2,'g--',tspan,ratio_real3,'m--') %Re(cyclic3) b/a
grid on;grid minor

ax4=subplot(2,3,4);
plot(ax4,tspan,ratio_im_cy1,'b-',tspan,ratio_im1,'r--',tspan,ratio_im2,'g--',tspan,ratio_im3,'m--') %Im(cyclic1) b/a
grid on;grid minor

ax5=subplot(2,3,5);
plot(ax5,tspan,ratio_im_cy2,'b-',tspan,ratio_im1,'r--',tspan,ratio_im2,'g--',tspan,ratio_im3,'m--') %Im(cyclic2) b/a
grid on;grid minor

ax6=subplot(2,3,6);
plot(ax6,tspan,ratio_im_cy3,'b-',tspan,ratio_im1,'r--',tspan,ratio_im2,'g--',tspan,ratio_im3,'m--') %Im(cyclic3) b/a
grid on;grid minor