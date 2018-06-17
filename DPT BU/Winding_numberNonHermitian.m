clc;clear;clf;
r=0.4;x=0;y=0.5;
N=1000;T=1;
N_k=400;
tspan=linspace(0,10*T,N);
time_evolved_states=zeros(N,2,N);
kcount=1;
options=odeset('RelTol',1e-13,'AbsTol',1e-14,'Refine',8); 
k_values=linspace(0,2*pi,N_k);

for k=linspace(0,2*pi,N_k)
%uses eigenstate as intial state %%%%%%%%%%%%%%%%%h2 or hBU%%%%%%%%%%%%%%%%%%%%%%%%
% [vec,val]=eig(hBU_matrix(0,r,k,0.5,0.7,T)); 
% initial_state=vec(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%h2 or hBU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uses cyclic state as initial state
[time1, state1]=ode45(@(t,y) H_BU(t,y,r,k,0,0.7,T),tspan,[1;0]);
[time2, state2]=ode45(@(t,y) H_BU(t,y,r,k,0,0.7,T),tspan,[0;1]);
u11=state1(end,1);
u12=state2(end,1);
u21=state1(end,2);
u22=state2(end,2);
[init_cyc,eigval] = eig([u11,u12;u21,u22]);
init_cyc1=init_cyc(:,1);
init_cyc2=init_cyc(:,2);
initial_state=init_cyc1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%h2 or hBU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[time,state]=ode45(@(t,psi) H_BU(t,psi,r,k,x,y,T),tspan,initial_state,options); 
time_evolved_states(:,:,kcount)=state;
kcount=kcount+1;
end

windingarray=zeros(1,N);
for time_index=1:N
time_stamp=tspan(time_index)
reduced_states_array=time_evolved_states(1:time_index,:,:);

geom_array=zeros(1,N);
geom_index=1;
for i=1:N_k %i denotes the microstate (i.e k value)
k_value=k_values(i);
ini_state=reduced_states_array(1,:,i);
LA=conj(ini_state)*transpose(reduced_states_array(end,:,i));
total_phase= real(-1i*log(LA/abs(LA)));

dyn_array=zeros(1,time_index);
dyn_index=1;
dt_prime=time_stamp/time_index;
    
    for j=1:time_index %j denotes the time values for a specified k
    t=time(j);
    final_state=reduced_states_array(j,:,i);
    %%%%%%%%%%%%%%%%%%%%%%h2 or hBU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dyn_array(dyn_index)=((conj(final_state)*hBU_matrix(t,r,k_value,x,y,T)*transpose(final_state))/(conj(final_state)*transpose(final_state)))*dt_prime; 
    dyn_index=dyn_index +1;
    end
    
dyn_phase=-real(sum(dyn_array));
geom_phase=total_phase-dyn_phase;
geom_array(geom_index)=geom_phase;
geom_index=geom_index+1;
end
% B=geom_array;
% B(isnan(geom_array))=0;
% B(isinf(geom_array))=0;
% % % % C=phasefix(B);
% % % % windingarray(time_index)=(C(end)-C(1))/(2*pi);
windingarray(time_index)=phasefix_windingcount(geom_array);
% subplot(1,2,1)
% plot(linspace(0,2*pi,N),geom_array,'b.')
% subplot(1,2,2)
% plot(linspace(0,2*pi,N),phasefix(B),'b.')
% plot(linspace(0,2*pi,N),phasefix(geom_array),'b.')
end
clf;
plot(tspan,windingarray,'linewidth',2)
% set(gca,'FontSize',35,'xtick',[2.05,2.92,4.97,9.73,9.93],'xticklabel',[],'ytick',[0,1,2,3,4,5])
% y=ylabel('\nu_{D}(t)','FontSize',50)
% x=xlabel('t','FontSize',50)
% set(x, 'Units', 'Normalized', 'Position', [0.5, -0.04, 0]);
% str={'2.05','2.92','4.97','9.73','9.93'}
% text([2.05-0.25,2.92-0.25,4.97-0.25,9.1,9.9],[-0.1,-0.1,-0.1,-0.1,-0.1],str,'fontsize',35)
set(gca,'FontSize',35,'xtick',[4.99,7.2,8.13,9.09],'xticklabel',[]);
y=ylabel('\nu_{D}(t)','FontSize',50);
x=xlabel('t','FontSize',50);
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.04, 0]);
str={'4.99','7.20','8.13','9.09'};
text([4.99-0.25,7.2-0.25,8.13-0.25,9.09-0.25],[-0.03,-0.03,-0.03,-0.03],str,'fontsize',35)