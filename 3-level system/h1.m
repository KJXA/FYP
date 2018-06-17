function hamiltonian = h1(t,mu,T)
hamiltonian = [0,1,0;1,0,0;0,0,0] + 1i*mu*(cos(2*pi*t/T)+1i)*[0,0,-1i;0,0,0;1i,0,0]+1i*mu*(sin(2*pi*t/T)+1i)*[0,0,0;0,0,1;0,1,0];
end