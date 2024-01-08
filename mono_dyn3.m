function [ dy ] = mono_dyn3(t,Y,k);
%system: LOTKA Volterra bi monomeric
% c1 <- monomer log(v)
% c2 <- monomer log(w)
% PARAMETERS
% k(1) = epsilon
% k(2) = c1

% Initialization
if ~iscolumn(Y)
    Yt = Y';
else
    Yt = Y;

c=Yt;


% ODE system
dc(1) = (k(1)-exp(c(2))-k(2)); %dc1/dt;
dc(2) = (exp(c(1))-k(1));
u=dc;

dy=u';


end