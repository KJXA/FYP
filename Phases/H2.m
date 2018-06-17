function psidot = H2(t,psi,ga,mu,T)
matrix = (ga*[1,0;0,-1] + 1i*mu*(cos(2*pi*t/T)+1i)*[0,1;1,0]);
psidot = (matrix * psi)/1i;
end