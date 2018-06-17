function psidot = H_2(t,psi,ga,mu,T)
matrix = (ga*[1,0,0;0,0,0;0,0,-1] + 1i*mu*(cos(2*pi*t/T)+1i)*[0,1,0;1,0,1;0,1,0]/sqrt(2))/1i;
psidot = matrix * psi;
end