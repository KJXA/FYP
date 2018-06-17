function final_H2=h2matrix_hop(t,r,k,x,y,T)
final_H2=(r*sin(k)+y)*[1,0;0,-1]+1j*(r*cos(k)+x)*(cos(2*pi*t/T)+1j)*[0,1;1,0];
end