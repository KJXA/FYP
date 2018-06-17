function matrix = h2_matrix(t,mu,T)
matrix = [1,0;0,-1] + 1i*mu*(cos(2*pi*t/T)+1i)*[0,1;1,0];
end