function psidot = H_BU(t,psi,r,k,x,y,T)
matrix = 1i*[0,1;(0.049)*exp(1i*2*pi*t/T)-(r*(cos(k)+1i*sin(k)) + x + 1i*y),0]; %set rho=0.3 or 0.049
psidot = (matrix * psi)/1i;
end
% 1i*[0,1;(0.049)*exp(1i*2*pi*t/T)-(sqrt(r)*(x+cos(k)+1i*(y+sin(k))))^2,0]
%old definition