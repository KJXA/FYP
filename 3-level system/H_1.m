function psidot = H_1(t,psi,mu,T)
matrix = ([0,1,0;1,0,0;0,0,0] + 1i*mu*(cos(2*pi*t/T)+1i)*[0,0,-1i;0,0,0;1i,0,0]+1i*mu*(sin(2*pi*t/T)+1i)*[0,0,0;0,0,1;0,1,0])/1i;
psidot = matrix * psi;
end