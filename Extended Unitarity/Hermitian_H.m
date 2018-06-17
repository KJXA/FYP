function psidot=Hermitian_H(t,psi,ga,mu,T)
matrix=[ga,mu*exp(-1i*2*pi*t/T);mu*exp(1i*2*pi*t/T),-ga];
psidot=(matrix*psi)/1i;
end