clc;clear;clf;
r=0.4;x=0;y=0.5;
N=1000;T=1; %T=1 for non-hopping; T=80 for hopping
tspan=linspace(0,10*T,N);
time_evolved_states=zeros(N,2,N);
kcount=1;
options=odeset('RelTol',1e-13,'AbsTol',1e-14,'Refine',8); %can lower to -10 -10

%model has period T=1; initial states = cyclic states of HBU at t=0, circle
%center at origin
for k=linspace(0,2*pi,N)
[t1,s1]=ode45(@(t,psi) H_BU(t,psi,r,k,0,0.7,T),[0,T],[1;0]);
[t2,s2]=ode45(@(t,psi) H_BU(t,psi,r,k,0,0.7,T),[0,T],[0;1]);
u11=s1(end,1);
u21=s1(end,2);
u12=s2(end,1);
u22=s2(end,2);
[cycvec,cycval]=eig([u11,u12;u21,u22]);
initial_state=cycvec(:,1);
[time,state]=ode45(@(t,psi) H_BU(t,psi,r,k,x,y,T),tspan,initial_state,options);
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
%     LE_array(LE_index) = dk*log((abs(LA)^2));
    LE_array(LE_index) = dk*log(((abs(LA)^2)/( sum(abs(time_evolved_states(i,:,j)).^2) )));
    LE_index = LE_index + 1;
    end
ratefn(ratefnindex) = sum(LE_array)/(-2*pi);
ratefnindex = ratefnindex + 1;
end
clf;
plot(tspan,ratefn,'b-','linewidth',2)
set(gca,'FontSize',35,'xtick',[4.99,7.2,8.13,9.09],'xticklabel',[]);
y=ylabel('g(t)','FontSize',50);
x=xlabel('t','FontSize',50);
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.04, 0]);
str={'4.99','7.20','8.13','9.09'};
text([4.99-0.25,7.2-0.25,8.13-0.25,9.09-0.25],[-1.2,-1.2,-1.2,-1.2],str,'fontsize',35)