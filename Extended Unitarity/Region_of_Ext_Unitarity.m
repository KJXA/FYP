clear;clc
%integration interval
tspan = linspace(0,1,100); 
%initial states
spin_up = [1;0];
spin_down = [0;1];
T=80;
options=odeset('RelTol',1e-13,'AbsTol',1e-14,'Refine',8);
tol = 0.05; %tolerance
N = 301; %no. of points
mu_array = zeros(1,N^2);
ga_array = zeros(1,N^2);
i = 1;
for ga= linspace(-15,15,N)
    for mu= linspace(-15,15,N)
        %get states at time t= 1
        [time1, state1] = ode45(@(t,psi) H1(t,psi,ga,mu),tspan,spin_up);
        [time2, state2] = ode45(@(t,psi) H1(t,psi,ga,mu),tspan,spin_down);
        
        %solve quadratic eigenvalue equation
        vec1 = state1(end,:); y11=vec1(1);y21=vec1(2);
        vec2 = state2(end,:); y12=vec2(1);y22=vec2(2);
        b = -y11-y22;
        c = y11*y22-y12*y21;
        eigen1 = (-b+sqrt(b^2-4*c))/2;
        eigen2 = (-b-sqrt(b^2-4*c))/2;
        
        %check if abs eigvals = 1
        if (abs(abs(eigen1)-1)< tol) && (abs(abs(eigen2)-1)< tol)
            mu_array(i) = mu;
            ga_array(i) = ga;
            i = i+ 1;
            
        end
    end
end
plot(mu_array, ga_array,'b.')
grid on
set(gca,'fontsize',25)
xlabel('\mu','FontSize',50)
ylabel('\gamma','FontSize',50)
axis([-15 15 -15 15])