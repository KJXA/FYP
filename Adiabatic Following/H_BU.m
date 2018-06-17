function psidot = H_BU(t,psi,rho,r,T)
matrix = 1i*[0,1;rho*exp(1i*2*pi*t/T)-r,0];
psidot = (matrix * psi)/1i;
end