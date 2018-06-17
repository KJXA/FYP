function psidot = H2_hop(t,psi,ga,mu,T)
matrix = (ga*[1,0;0,-1] + 1i*mu*(cos(2*pi*t/T)+1i)*[0,1;1,0])/1i;
psidot = matrix * psi;
