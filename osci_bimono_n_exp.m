function [ dy ] = osci_bimono_n_exp(t,Y,num_poly);
%system: y(1)=mono, y(2)=mono2, y(3:end)=polymers
%tc1+c(k)->c(k+1) and c(k+1)+c1 -> c(k) + 2c1
%Y(1)=c1, Y(2)=tc1, Y(3)=c2, Y(4)=c3
%num_poly: size of the vector y - 2 = number of sizes taken by the polymers
%there are (2*num_poly-1) reactions, num_poly+2 species, 4*num_poly-2 complexes


xx = Y(1);
xy = Y(2);
c(1:num_poly) = Y(3:num_poly+2);
eps = sum(c(1:num_poly));
u(1) = eps - c(1) -exp(xy); %da1/dt
u(2) = exp(xx) - (eps-c(num_poly)); %db1/dt
dc(1) = -c(1)*exp(xy)+c(2)*exp(xx); %dc1/dt;
dc(2:num_poly-1) = c(3:num_poly).*exp(xx)...
                   + c(1:num_poly-2).*exp(xy)...
                   - c(2:num_poly-1).*(exp(xx)+exp(xy)); %dci/dt
dc(num_poly) = c(num_poly-1)*exp(xy)-c(num_poly)*exp(xx); %dcn/dt
u(3:num_poly+2) = dc(1:num_poly);

dy = u';


end
