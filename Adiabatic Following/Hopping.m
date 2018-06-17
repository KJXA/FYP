clear;clc;clf
%Number of points to evaluate
N = 10000;
ga=1; %rho for BU model
T=80;
mu=0.2; %hopping starts at around 0.81(H2) %r for BU model; hopping starts at rho/r=0.43 (BU)
%integration interval
tspan = linspace(0,T,N); %4 periods still ok; 5 periods not ok

%initial states
spin_up = [1;0];
spin_down = [0;1];
options=odeset('RelTol',1e-16,'AbsTol',1e-25,'Refine',8);
%get cyclic states
[time1, state1]=ode45(@(t,y) H2(t,y,ga,mu,T),tspan,spin_up,options);
[time2, state2]=ode45(@(t,y) H2(t,y,ga,mu,T),tspan,spin_down,options);

cyclic_states1 = zeros(2,N);
cyclic_states2 = zeros(2,N);
u11=state1(end,1);
u12=state2(end,1);
u21=state1(end,2);
u22=state2(end,2);
[init_cyc,eigval] = eig([u11,u12;u21,u22]);
init_cyc1=init_cyc(:,1);
init_cyc2=init_cyc(:,2);
for i = 1:N
    u11=state1(i,1);
    u12=state2(i,1);
    u21=state1(i,2);
    u22=state2(i,2);
    cyclic_states1(:,i)=[u11,u12;u21,u22]*init_cyc1;
    cyclic_states2(:,i)=[u11,u12;u21,u22]*init_cyc2;
end

ratio_real_cy1 = real(cyclic_states1(2,:)./cyclic_states1(1,:));
ratio_im_cy1 = imag(cyclic_states1(2,:)./cyclic_states1(1,:));
ratio_real_cy2 = real(cyclic_states2(2,:)./cyclic_states2(1,:));
ratio_im_cy2 = imag(cyclic_states2(2,:)./cyclic_states2(1,:));

%get instantaneous eigenstates
invec1 = zeros(2,N);
invec2 = zeros(2,N);
index = 1;
for t = tspan
[invec,invalue] = eig(h2_matrix(t,mu,T)); %rho=ga;r=mu;
invec1(:,index)=invec(:,1);
invec2(:,index)=invec(:,2);
index = index + 1;
end

ratio_real1 = real(invec1(2,:)./invec1(1,:));
ratio_im1 = imag(invec1(2,:)./invec1(1,:));
ratio_real2 = real(invec2(2,:)./invec2(1,:));
ratio_im2 = imag(invec2(2,:)./invec2(1,:));
%plot ratios

ax1=subplot(1,2,1);
plot(ax1,tspan,ratio_real_cy1,'b-',tspan,ratio_real1,'r--',tspan,ratio_real2,'r--','linewidth',2) %Re(cyclic1)
grid on;grid minor
set(gca,'FontSize',30)
xlabel('t','FontSize',30)
ylabel('Re(\psi_+)','FontSize',30)
% ax2=subplot(1,2,1);
% plot(ax2,tspan,ratio_real_cy2,'b-',tspan,ratio_real1,'r--',tspan,ratio_real2,'r--','linewidth',2) %Re(cyclic2)
% grid on;grid minor
% set(gca,'FontSize',30)
% xlabel('t','FontSize',30)
% ylabel('Re(\psi_-)','FontSize',30)
ax3=subplot(1,2,2);
plot(ax3,tspan,ratio_im_cy1,'b-',tspan,ratio_im1,'r--',tspan,ratio_im2,'r--','linewidth',2) %Im(cyclic1)
grid on;grid minor
set(gca,'FontSize',30)
xlabel('t','FontSize',30)
ylabel('Im(\psi_+)','FontSize',30)
% ax4=subplot(1,2,2);
% plot(ax4,tspan,ratio_im_cy2,'b-',tspan,ratio_im1,'r--',tspan,ratio_im2,'r--','linewidth',2) %Im(cyclic2)
% grid on;grid minor
% set(gca,'FontSize',30)
% xlabel('t','FontSize',30)
% ylabel('Im(\psi_-)','FontSize',30)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot for thesis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(tspan,ratio_im_cy1,'b-',tspan,ratio_im1,'r--',tspan,ratio_im2,'r--','linewidth',2) 
% plot(tspan,ratio_im_cy2,'b-',tspan,ratio_im1,'r--',tspan,ratio_im2,'r--','linewidth',2) 
% grid on;grid minor
% set(gca,'FontSize',50)
% xlabel('t','FontSize',50)
% ylabel('Im(\psi_-)','FontSize',50)