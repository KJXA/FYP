function psidot = H1(t,psi,ga,mu)
matrix = (ga*[1,0;0,-1]+mu*1i*(cos(2*pi*t)+sin(4*pi*t))*[0,1;1,0])/(1i);
psidot = matrix*psi;
end