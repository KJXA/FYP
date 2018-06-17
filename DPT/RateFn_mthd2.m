clc;clear;clf;
mu=0.6;r=0.3;ga=1;
N=1000;
tspan=linspace(0,20,1000);
dk=2*pi/N;
ratefnindex=1;
rate_fn=zeros(1,N);
for t= linspace(0,20,1000)
    LE_array=zeros(1,N);
    LEindex=1;
    for k = linspace(0,2*pi,N)
    [vec,values] = eig(H_f(mu,r,k,0));
    initial_state1 = vec(:,1);
    initial_state2 = vec(:,2);
    LA_1 = transpose(conj(initial_state1))*expm(-1i*t*H_f(mu,r,k,ga))*initial_state1;
%     LE_array(LEindex) = dk*log((abs(LA_1)^2)/sum(abs(expm(-1i*t*H_f(mu,r,k,ga))*initial_state1).^2));
    LE_array(LEindex) = dk*log((abs(LA_1)^2)/abs(conj(transpose(expm(-1i*t*H_f(mu,r,k,ga))*initial_state1))*(expm(-1i*t*H_f(mu,r,k,ga))*initial_state1)));
    LEindex = LEindex + 1;
    end
    rate_fn(ratefnindex) = -1/(2*pi)*sum(LE_array);
    ratefnindex = ratefnindex + 1; 
end
clf;
plot(tspan,rate_fn,'b-','linewidth',2)
set(gca,'FontSize',50,'xtick',[4.7,8.7,12.8,16.9],'xticklabel',{4.7,8.7,12.8,16.9})
xlabel('t','FontSize',50)
ylabel('g(t)','FontSize',50)
ylim([0,1.2])