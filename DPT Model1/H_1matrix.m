function matrix = H_1matrix(t,r,k,x,y,T)
matrix = ((r*sin(k)+y)*[1,0;0,-1]+(r*cos(k)+x)*1i*cos(2*pi*t)*[0,1;1,0]);
end