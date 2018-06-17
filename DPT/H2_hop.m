function psidot=H2_hop(t,psi,r,k,x,y,T)
matrix = (r*sin(k)+y)*[1,0;0,-1]+1j*(r*cos(k)+x)*(cos(2*pi*t/T)+1j)*[0,1;1,0];
psidot = (matrix*psi)/1j;
end