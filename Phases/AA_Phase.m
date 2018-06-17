clear;clc;%clf
%Number of points to evaluate
N = 10000;
ga=1;
T=50;
mu=1.2;
%integration interval
tspan = linspace(0,T,N); 

% initial states
spin_up = [1;0];
spin_down = [0;1];
options=odeset('RelTol',1e-16,'AbsTol',1e-25,'Refine',16);
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

% % hamiltonian at each point in time
% hamiltonians=zeros(2,2,10);
% for i=1:N
%     hamiltonians(:,:,i)=h2_matrix(tspan(i),mu,T);
% end

dynamical_array=zeros(N,1);
for i=1:N
dynamical_array(i)= (transpose(conj(cyclic_states1(:,i)))*h2_matrix(tspan(i),mu,T)*cyclic_states1(:,i))/(sum(abs(cyclic_states1(:,i)).^2));
end
negative_dynamical_phase=real(sum(dynamical_array))

alpha=(-1j*log(cyclic_states1(1,end)/cyclic_states1(1,1)))
AA_phase=mod(real(alpha)+negative_dynamical_phase,2*pi)


