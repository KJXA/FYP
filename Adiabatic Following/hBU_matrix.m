function matrix = hBU_matrix(t,rho,r,T)
matrix = 1i*[0,1;rho*exp(1i*2*pi*t/T)-r,0];
end