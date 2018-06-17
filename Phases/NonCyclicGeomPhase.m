clear;clc;clf
%Number of points to evaluate
N = 200;
ga=1;
T=80;
n=100;
geom_phase_array=zeros(1,n);

mu=1.0;
% Possible initial states to try
spin_up = [1;0];
spin_down = [0;1];
[vecs,vals]=eig(h2_matrix(0,mu,T));
eig1=vecs(:,1);
eig2=vecs(:,2);
%Choose initial state
initial_state=eig1;
%T=80; choose time parameter to evolve state till
geom_array=zeros(1,N);
geom_index=1;

for t=linspace(0.001,T,N);
% t=40; %end point
tspan=linspace(0,t,N);
%get time-evolved states from t=0 to t=t
options=odeset('RelTol',1e-16,'AbsTol',1e-25,'Refine',8);
[time,state]=ode45(@(t,y) H2(t,y,ga,mu,T),tspan,initial_state,options);
%get dynamical phase
dyn_array=zeros(1,length(state));
dyn_index=1;
dt_prime=t/length(state);
for t_prime=tspan
    final_state=state(dyn_index,:);
    d_dyn=conj(final_state)*h2_matrix(t_prime,mu,T)*transpose(final_state)*dt_prime/(conj(final_state)*transpose(final_state));
    dyn_array(dyn_index)=d_dyn;
    dyn_index=dyn_index+1;
end
dyn_phase=-real(sum(dyn_array));
%get Loschmidt amplitude and thus total phase
LA=transpose(conj(initial_state))*transpose(final_state);
total_phase=real(-1i*log(LA/abs(LA))); %total phase real; might have numerical errors
geom_phase=total_phase-dyn_phase;
geom_array(geom_index)=geom_phase;
geom_index=geom_index+1;
end

new_geom_array=phasefix(geom_array);
subplot(1,2,1)
plot(tspan,geom_array,'b.')
subplot(1,2,2)
plot(tspan,phasefix(phasefix(new_geom_array)),'b.')