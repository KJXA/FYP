function psidot = H2(t,psi,ga,mu)
matrix = (ga*[1,0;0,-1]+1i*mu*(sin(2*pi*t)+1i)*[0,1;1,0])/1i;
psidot = matrix*psi;
end