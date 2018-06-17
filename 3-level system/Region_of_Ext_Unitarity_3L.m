clear;clc
% initial states
spin_up = [1;0;0];
spin_mid = [0;1;0];
spin_down = [0;0;1];

tol = 0.05; %tolerance
N = 101; %no. of points
T=100;

%integration interval
tspan = linspace(0,1,T);

mu_array = zeros(1,N);
ga_array = zeros(1,N);
i = 1;
for ga= linspace(-15,15,N)
    for mu= linspace(-15,15,N)
        [time1, state1] = ode45(@(t,psi) H2ga_mu(t,psi,ga,mu,T),tspan,spin_up);
        [time2, state2] = ode45(@(t,psi) H2ga_mu(t,psi,ga,mu,T),tspan,spin_mid);
        [time3, state3] = ode45(@(t,psi) H2ga_mu(t,psi,ga,mu,T),tspan,spin_down);
        
        % get eigenvalues
        temp1=transpose(state1(end,:));
        temp2=transpose(state2(end,:));
        temp3=transpose(state3(end,:));
        floquet = [temp1 temp2 temp3];
        [vec,eigen] = eig(floquet,'vector');
        eigen = abs(eigen).^2;
        %check if abs eigvals = 1
        if (abs(eigen(1)-1)< tol) && (abs(eigen(2)-1)< tol) && (abs(eigen(3)-1)< tol)
            mu_array(i) = mu;
            ga_array(i) = ga;
            i = i+ 1;
            
        end
    end
end
plot(mu_array, ga_array,'b.')
grid on
xlabel('\mu')
ylabel('\gamma')
axis([-15 15 -15 15])