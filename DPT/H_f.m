function final_H = H_f(mu,r,k,ga)
final_H=(mu+r*cos(k))*[0,1;1,0]+(r*sin(k)+1j*ga/2)*[1,0;0,-1];
end