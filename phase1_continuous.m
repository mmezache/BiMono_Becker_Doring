% Illustration of the bimonomeric pol/depol model
% phase 1 Continuous approximation

clear all
close all


%% Setting parameters
L = 250; % Maximum size of clusters
dx = 0.5; % size step
sig = 10; % spread of the size distribution
T = 10000; % final time
dt = 0.05; % time step
NX = floor (L/dx ) - 1;
xsize = linspace( 0, L, floor( L/dx ) );
size_distr = zeros( floor(T/dt) ,500 );
size_distr( 1, : ) = half_gaussian( xsize, sig, 0.02 );

eps = dx*sum( size_distr(1,:) ); % total concentration
time_span = linspace( 0, T, floor( T/dt ) );
k = [eps 0]; % parameters for the LV system; 
%% Compute monomers
Y0  = zeros( 2, 1 );
Y0(1) = log(0.6);
Y0(2) = log(0.6);

[Tout, Yout] = ode89( @mono_dyn3, time_span, Y0, [], k );

v = exp( Yout( :, 1 ) );
w = exp( Yout(:, 2 ) );

%% Compute PDE
% Implicit scheme 
% Dirichlet Neumann boundary conditions
dd = v + w;
vv = w - v;


a = -( dt/( 2*dx ) )*vv - ( dt/( 2*dx*dx ) )*dd;
b = 1 + dd*( dt/( dx*dx ) );
c = ( dt/( 2*dx ) )*vv - ( dt/( 2*dx*dx ) )*dd;

r =dd./( 2*dx*vv + dd );
resol = zeros( NX ,1 );
for ie =2:floor( T/dt )
    Q = diag( b(ie)*ones( 1, NX ) ) + diag( c(ie)*ones( 1, NX-1 ), 1 ) + diag( a(ie)*ones( 1, NX-1 ), -1 );
    Q(1,1) = 1 - a(ie);
    resol = size_distr( ie-1, 1:NX );
    size_distr( ie, 1:NX ) = tridiagonal( Q, resol' )';

end
%% Conserved Quantities 

eps_toto = sum( size_distr, 2 )*dx;
mtot = v + w + dx*( xsize*size_distr' )';

figure 
subplot(1, 2, 1)
plot( Tout, eps_toto )
title( 'Total number of polymers' )
subplot(1, 2, 2)
plot( Tout, mtot )
title( 'Total mass' )

% explicitation for the total mass
temp_vw = ( v(2:end) - v(1:end-1) )./dt + ( w(2:end) - w(1:end-1) )./dt;
temp_dist = ( v(2:end) - w(2:end) ).*( eps_toto(2:end) );
temp_c0 = ( ( 1 + dx )*v + ( 1 - dx )*w ).*( size_distr( :, 1 )./2 );
figure 
subplot(1, 2, 1)
plot( Tout(2:end), temp_vw - temp_dist )
title( 'Diffenrence between residuals without the contribution of c_0' )
subplot(1, 2, 2)
plot( Tout, temp_c0 )
title( 'Contribution of the bound on the mass' )

%% Characteristic quantities
n_iter = length(eps_toto);
ENER = zeros( n_iter, 1 );

for ie = 1:n_iter
    ENER(ie) = v(ie) + w(ie) -2*( eps_toto(ie) ) ...
               - eps_toto(ie)*log( ( v(ie)*w(ie) )/( eps_toto(ie) )^2 );
end


%% Animation size distribution
N = 500; % length of the vector
nframes = 400; % number of frames

figh = figure( 'position', [100 100 850 600] );
for k = 1:nframes
    % Wipe the slate clean so we are plotting with a black figure
    clf
    % Extract Data 
    v_k = v( 1 + (k-1)*N );
    w_k = w( 1 + (k-1)*N );
    y_k = size_distr( 1 + (k-1)*N, : );
    
    subplot(1,2,2)
    plot( xsize,y_k, 'LineWidth', 2 );
    axis( [0 200 0 max( max( size_distr(1, :) ) , max( size_distr( end, :) ) ) ] )
    grid on
    
    
    
    subplot(1, 2, 1)
    % Plot the current location of the particle
    plot( v_k, w_k, 'ro', 'LineWidth', 3, 'MarkerSize', 15 )

    % Plot the entire trajecory
    hold on
    plot(v, w, 'b')
    
    % Decorate the plot
    grid on;
    xlabel( 'V' );
    ylabel( 'W' );
    axis( [0 max(v) 0 max(w)] )

    % Force Matlab to draw the image at this point
    movFrame(k) = getframe(figh);
    
end

%% Figure about the change of the shape of the size distribution
index_shape = [27324 81970 136615 191261];
figure1 = figure;
cmap = copper( length(index_shape) );
grid on
for ie = 1:length(index_shape)
    plot( xsize(1:400), size_distr( index_shape(ie), 1:400 ) )
    hold on
end
hcb = colorbar;

set(gca, 'colororder', cmap, 'colormap', cmap)
caxis( [Tout( index_shape(1) ) Tout( index_shape(end) ) ] )
% set ticklabel length
tl = split( sprintf('%.0f\n',hcb.Ticks) );
hcb.TickLabels = tl(1:end-1);
% set cb title
title( hcb, 'Time' )
xlabel( 'size' )

index_shape2 = [1 54647 109293 163938];
figure2 = figure;
cmap = copper( length(index_shape2) );
grid on
for ie = 1:length(index_shape2)
    plot( xsize(1:100), size_distr( index_shape2(ie), 1:100 ) )
    hold on
end
hcb = colorbar;
set(gca, 'colororder', cmap, 'colormap' ,cmap)
caxis( [Tout( index_shape2(1) ) Tout( index_shape2(end) ) ] )
% set ticklabel length
tl = split( sprintf('%.0f\n', hcb.Ticks) );
hcb.TickLabels = tl(1:end-1);
% set cb title
title( hcb, 'Time' )
xlabel( 'size' )

%% Additional Functions 
function [y] = half_gaussian(x, sig, m)
    y = exp( - ( (x).^2 )/(2*sig) )*sqrt( 2/( pi*sig ) )*m*(1.1262/1.1262);
end
