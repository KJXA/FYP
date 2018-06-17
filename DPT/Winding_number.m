clc;clear;clf;
mu=0.6;r=0.3;ga=1;
N=300;
tspan=linspace(0,20,N);
dk=2*pi/N;
windingarray=zeros(1,N);
for i=1:N
t=tspan(i);
% k=pi
geom_array=zeros(1,N);
geom_index=1;
for k = linspace(0,2*pi,N)
    [vec,values] = eig(H_f(mu,r,k,0));
    initial_state1 = vec(:,1);
    initial_state2 = vec(:,2);
    LA_1 = transpose(conj(initial_state1))*expm(-1i*t*H_f(mu,r,k,ga))*initial_state1;
    total_phase= real(-1i*log(LA_1/abs(LA_1))); %total phase is real but imag part non-zero numerically
    
    dyn_array=zeros(1,N);
    dyn_index=1;
    dt_prime=t/N;
    for t_prime=linspace(0,t,N)
        final_state1=expm(-1i*t_prime*H_f(mu,r,k,ga))*initial_state1; %final state at time=t'
        final_state2=expm(-1i*t_prime*H_f(mu,r,k,ga))*initial_state2;
        dyn_1=(transpose(conj(final_state1))*H_f(mu,r,k,ga)*final_state1)/(transpose(conj(final_state1))*final_state1)*dt_prime;
        dyn_array(dyn_index)=dyn_1;
        dyn_index=dyn_index +1;
    end
%     definition 1 of dyn phase
    dyn_phase=-real(sum(dyn_array));
%     definition 2 of dyn phase
%     final1=expm(-1i*t*H_f(mu,r,k,ga))*initial_state1; %final state at time=t
%     dyn_phase2=-sum(dyn_array)+1i*log(sqrt((transpose(conj(final1))*final1)/(transpose(conj(initial_state1))*initial_state1)));
%     definition 1 used
    geom_phase=total_phase-dyn_phase;
    geom_array(geom_index)=geom_phase;
    geom_index=geom_index+1;
end
temp=phasefix(geom_array);
winding=(temp(end)-temp(1))/(2*pi);
windingarray(i)=winding;
end
% subplot(1,2,1)
% plot(linspace(0,2*pi,N),geom_array,'b.')
% subplot(1,2,2)
% plot(linspace(0,2*pi,N),phasefix(geom_array),'b-')

clf;
plot(tspan,windingarray,'b-','linewidth',2)
set(gca,'FontSize',50,'ytick',[0 1 2 3 4],'xtick',[4.7,8.7,12.8,17.2],'xticklabel',{4.7,8.7,12.8,16.9})
xlabel('t','FontSize',50)
ylabel('\nu_{D}(t)','FontSize',50)
