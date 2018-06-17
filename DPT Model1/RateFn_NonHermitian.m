clc;clear;clf;
r=0.5;x=2.5;y=2.5;
N=1000;T=1; %T=1 for non-hopping; T=80 for hopping
tspan=linspace(0,10*T,N);
time_evolved_states=zeros(N,2,N);
kcount=1;
options=odeset('RelTol',1e-14,'AbsTol',1e-14,'Refine',8); %can lower to -10 -10
%initial states = eigvecs of H2 at t=0, circle center at origin 
% for k=linspace(0,2*pi,N)
% [vec,val]=eig(H_1matrix(0,r,k,5,0.5,T));
% initial_state=vec(:,1);
% [time,state]=ode45(@(t,psi) H_1(t,psi,r,k,x,y,T),tspan,initial_state); 
% time_evolved_states(:,:,kcount)=state;
% kcount=kcount+1;
% end

%model has period T=1; initial states = cyclic states of H2 at t=0, circle
%center at origin
for k=linspace(0,2*pi,N)
[t1,s1]=ode45(@(t,psi) H_1(t,psi,r,k,0,0,T),[0,T],[1;0]);
[t2,s2]=ode45(@(t,psi) H_1(t,psi,r,k,0,0,T),[0,T],[0;1]);
u11=s1(end,1);
u21=s1(end,2);
u12=s2(end,1);
u22=s2(end,2);
[cycvec,cycval]=eig([u11,u12;u21,u22]);
initial_state=cycvec(:,1);
[time,state]=ode45(@(t,psi) H_1(t,psi,r,k,x,y,T),tspan,initial_state);
time_evolved_states(:,:,kcount)=state;
kcount=kcount+1;
end

dk=2*pi/N;
ratefn=zeros(1,N);
ratefnindex=1;

for i=1:N %iterate over time
LE_array = zeros(1,N);
LE_index = 1;
    for j=1:N %iterate over k
    initial_psi=time_evolved_states(1,:,j);
    LA = conj(initial_psi)*transpose(time_evolved_states(i,:,j));
    LE_array(LE_index) = dk*log(abs(LA)^2);
    LE_index = LE_index + 1;
    end
ratefn(ratefnindex) = sum(LE_array)/(-2*pi);
ratefnindex = ratefnindex + 1;
end

plot(tspan,ratefn,'b-')