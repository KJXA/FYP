function psidot = H_BU(t,psi,r,k,x,y,T)
matrix = 1i*[0,1;(r*sin(k)+y)*exp(1i*2*pi*t/T)-(r*cos(k)+x),0];
psidot = (matrix * psi)/1i;
end