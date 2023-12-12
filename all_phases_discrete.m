% Illustration of the bimonomeric pol/depol model 
% Complete numerical simulations of the discrete system

clear all
close all

%% Setting parameters
num_poly = 500; %number of sizes of polymers
Y0 = zeros(num_poly+2,1); % Allocate and name the initial conditions

% Monomers 
Y0(1) = log(0.6);
Y0(2) = log(0.6);

% Parameters of the initial distribution 
sig = 10; % spread
eps = 0.02; % sum of the cluster distribution
xsize = linspace( 1, num_poly, num_poly ); % size range

% Initial distribution 
Y0(3:num_poly+2) = half_gaussian(xsize,sig, eps)*1.1493;

%% ODE Computations
T0 = 0; Tf = 300000;
dt = 0.5;

[Tout,Yout] = ode89( @osci_bimono_n_exp, [T0 Tf], Y0, [], num_poly );
v = exp(Yout(:,1));
w = exp(Yout(:,2));

%% Concerved quantities

masstot = v + w + ( (1:num_poly)*Yout( :, 3:num_poly+2 )' )';
num_tot = sum( Yout(:,3:end)' );

% Plot figures
figure    
subplot(211)
plot(Tout,masstot)
title(['total mass (should be constant), num polymers=',num2str(num_poly)])
subplot(212)
plot(Tout,num_tot);
title('epsilon (should be constant)')

%% Characteristic quantities

n_iter = length(Yout(:,1));
LEN = zeros(n_iter,1);
SPREAD = zeros(n_iter,1);
ENER = zeros(n_iter,1);

% Length, Spread, Energy
for ie = 1:n_iter 
    LEN(ie) = ( (1:num_poly)*Yout(ie,3:num_poly+2)' )/sum( Yout(ie,3:num_poly+2) );
    SPREAD(ie) = ( ( ( (1:num_poly) - LEN(ie) ).^2 )*Yout(ie,3:num_poly+2)' )/sum( Yout(ie,3:num_poly+2) );
    ENER(ie) = v(ie) + w(ie) -2*( ones(1,num_poly)*Yout(ie,3:num_poly+2)' ) ...
               - ( ones(1,num_poly)*Yout(ie,3:num_poly+2)' )...
               *log( ( v(ie)*w(ie) )/( ones(1,num_poly)*Yout(ie,3:num_poly+2)' )^2 );
end


figure

% Top two plots
tiledlayout(2,2)
nexttile
plot(Tout,LEN)
title('Evolution of the characteristic length')
nexttile
plot(Tout,SPREAD)
title('Evolution of the characteristic spread')
% Plot that spans
nexttile([1 2])
semilogy(Tout,ENER)
hold on 
yline(eps,Color='red')
yline(eps*eps,Color='red')
hold off
title('Evolution of the characteristic energy')

%% Poincare section
n_iter = 22457;
temp_v = v(5337:27793);
temp_w = w(5337:27793);
L1 = [temp_v';temp_w'];
temp = linspace(0.7,0.02,n_iter);
L2 = [temp;temp-Yout(end,3)];
P_inter = InterX(L1,L2); % sequence of R^2 = points where trajectory cut the poincare map
t_trun = Tout(5337:27793);

%% Plot Phase Space
figure;
cmap = jet(130);
plot( temp_v, temp_w, "k" )
for ie = 1:129
    hold on
    plot( P_inter( 1, 130-ie+1 ), P_inter( 2, 130-ie+1 ), 'o', 'LineWidth', 3, 'MarkerSize', 15 )
end
hcb = colorbar;
set(gca, 'colororder', cmap, 'colormap', cmap)
caxis([t_trun(1) t_trun(end)])
% set ticklabel length
tl = split( sprintf('%.0f\n', hcb.Ticks) );
hcb.TickLabels = tl(1:end-1);
% set cb title
title(hcb, 'Time')
xlabel('v')
ylabel('w')

%% Plot Energy decay
x_ener = linspace(5337, 27793, 130);
figure;
cmap = jet(130);
for ie = 1:130
    plot( Tout( floor( x_ener(ie) ) ), ENER( floor( x_ener(ie) ) ), 'o', 'LineWidth', 3, 'MarkerSize', 15 )
    hold on
end
hcb = colorbar;
set(gca, 'colororder', cmap, 'colormap', cmap)
caxis( [t_trun(1) t_trun(end)] )
% set ticklabel length
tl = split( sprintf( '%.0f\n', hcb.Ticks ) );
hcb.TickLabels = tl(1:end-1);
% set cb title
title( hcb, 'Time' )
plot( Tout, ENER, "k" )
xlabel('time')
ylabel('Energy')

%% Animation size distribution

N = length( Yout(:,1) ); % length of the vector
tin = 1; % initial time index
tout = 5337; % final time index
nframes = tout - tin; % number of frames
toto2 = sum( Yout(1,3:end) ); % epsilon (illustrate the phase space)
n_gap = 5; % number of frames that are not plotted
figh = figure( 'position', [100 100 850 600] );
for k = 1:nframes
    % Wipe the slate clean so we are plotting with a black figure
    clf
    % Extract Data 
    v_k = v( n_gap*(tin +k) );
    w_k = w( n_gap*(tin +k) );
    y_k = Yout( n_gap*(tin +k), 3:num_poly+2 );
    En_k = real( ENER( n_gap*(tin +k) ) );
    Sig_k = SPREAD( n_gap*( tin +k ) );
    T_k = Tout( n_gap*( tin +k ) );
    

    subplot(2, 3, [4 5 6])
    plot(y_k, 'LineWidth', 2);
    axis([0 num_poly 0 max( max( Yout(:,3:num_poly+2) ) )])
    yline(toto2^2,'r-')
    grid on
    
    subplot(2,3,1)
    % Plot the current location of the particle
    plot(v_k, w_k, 'ro', 'LineWidth', 3, 'MarkerSize', 15)

    % Plot the entire trajecory
    hold on
    plot( v( n_gap*tin:n_gap*tout ), w( n_gap*tin:n_gap*tout ), 'b')
    
    % Decorate the plot
    grid on;
    xlabel('V');
    ylabel('W');
    xline(toto2, 'r-')
    yline(toto2 - Yout(end,3), 'r-')
    axis([0 max(v) 0 max(w)])

    subplot(2,3,2)
    % Plot the current location of the particle
    plot(T_k, En_k, 'ro', 'LineWidth', 3, 'MarkerSize', 15)

    % Plot the entire trajecory
    hold on
    semilogy( Tout( n_gap*tin:n_gap*tout-2 ), real( ENER( n_gap*tin:n_gap*tout-2 ) ), 'b')
    
    % Decorate the plot
    grid on;
    xlabel('Time');
    ylabel('Energy');
    %axis([0 max(Tout) 0 max(ENER)])
    
    subplot(2,3,3)
    % Plot the current location of the particle
    plot(T_k, Sig_k, 'ro', 'LineWidth', 3, 'MarkerSize', 15)

    % Plot the entire trajecory
    hold on
    plot(Tout( n_gap*tin:n_gap*tout-2 ),SPREAD( n_gap*tin:n_gap*tout-2 ), 'b')
    
    % Decorate the plot
    grid on;
    xlabel('Time');
    ylabel('Spread');
    


    % Force Matlab to draw the image at this point
    %drawnow
    movFrame(k) = getframe(figh);
    
end


%% Additional Functions 
function [y]=half_gaussian(x,sig,m)
    y = exp(-((x).^2)/(2*sig))*sqrt(2/(pi*sig))*(m)*(1.1262/1.1262);
end
function [y]=gaussian(x,m,sig,eps);
    y = exp(-((x-m).^2)/(2*sig))*(1/sqrt(2*pi*sig))*eps ;
end
