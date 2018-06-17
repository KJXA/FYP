function matrix = hBU_matrix(t,r,k,x,y,T)
matrix = 1i*[0,1;(0.049)*exp(1i*2*pi*t/T)-(r*(cos(k)+1i*sin(k)) + x + 1i*y),0]; %set rho=0.3 or 0.049
end