function psidot = H_test(t,psi,ga,mu)
matrix = (ga*[1,0;0,-1]+mu*1i*sin(t)*[0,1;1,0])/1i;
psidot = matrix*psi;
end