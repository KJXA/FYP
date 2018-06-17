clear;clc;clf
%Number of points to evaluate
N = 2000;
ga=1;
T=80;
n=100;
geom_phase_array=zeros(1,n);
options=odeset('RelTol',1e-16,'AbsTol',1e-25,'Refine',8);
mu=1.0;
% Possible initial states to try
spin_up = [1;0];
spin_down = [0;1];
[vecs,vals]=eig(h2_matrix(0,mu,T));
eig1=vecs(:,1);
eig2=vecs(:,2);
%%%%%%%%%%%%
[time1,s1]=ode45(@(t,y) H2(t,y,ga,mu,T),[0,T],[1;0],options);
[time2,s2]=ode45(@(t,y) H2(t,y,ga,mu,T),[0,T],[0;1],options);
u11=s1(end,1);
u21=s1(end,2);
u12=s2(end,1);
u22=s2(end,2);
[cycvec,cycval]=eig([u11,u12;u21,u22]);
cy1=cycvec(:,1);
cy2=cycvec(:,2);
%Choose initial state
initial_state=cy1;
T=80; %choose time parameter to evolve state till
geom_array=zeros(1,N);
geom_index=1;

tspan=linspace(0,T,N);
%get time-evolved states from t=0 to t=T
[time,state]=ode45(@(t,y) H2(t,y,ga,mu,T),tspan,initial_state,options);
%get dynamical phase
% dyn_array=zeros(1,length(state));
% dyn_index=1;
dt=T/length(time);
for i=linspace(1,length(time),N)
    final_state=state(i,:);
    LA=transpose(conj(initial_state))*transpose(final_state);
    total_phase=real(-1i*log(LA/abs(LA))); %total phase real;
    dyn_array=zeros(1,i);
    for dyn_index=1:i
        evol_state=state(dyn_index,:);
        dyn_array(dyn_index)=conj(evol_state)*h2_matrix(time(dyn_index),mu,T)*transpose(evol_state)/(conj(evol_state)*transpose(evol_state))*dt;
    end
    geom_phase=total_phase+real(sum(dyn_array));
    geom_array(geom_index)=geom_phase;
    geom_index=geom_index +1;
end

new_geom_array=phasefix(geom_array);
subplot(1,2,1)
plot(tspan,geom_array,'b.')
subplot(1,2,2)
plot(tspan,new_geom_array,'b.')
% new_geom_array(99:end)=new_geom_array(99:end)+2;
