%% |Spatio-Temporal Calculation|
%% |Adapted from read_data.m, generated form Kundu et. al, written by Dr. Tamas Jozsa|
% |Start|

clear
close all
clc
% |READ CHANNEL PARAMETERS AND SPATIO-TEMPORAL AVERAGES|

% read parameters
% Lz, Lx, Ly, nu, Delta p
params = xlsread('Reynolds_stresses.xlsx','parameters');

Lz = params(1); Lx = params(2); Ly = params(3);
nu = params(4); % kinematic viscosity

% these values are equal to unity because they are the reference quantities
% used to make the data dimensionless
u_b = 1.0; % bulk velocity (average velocity in the entire channel)
rho = 1.0; % density
delta = Lx/2; % boundary layer thickness  =  channel half-height

% bulk Reynolds number based on channel half height and mean velocity
Re_b  =  u_b*delta/nu;

% read wall-normal coordinate and spatio-temporal averages
% x, <w>, <w'w'>, <u'u'> , <v'v'>, <u'w'>
ST_ave_dat = xlsread('Reynolds_stresses.xlsx','Reynolds_stresses');


%% READ TIME SAMPLES AT PROBES PLACED ALONG A WALL_NORMAL LINE

hinfo  =  hdf5info('time_samples.hdf5');

% sampling time
t_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(1));
% wall-normal location of the samples
x_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(5))+1.0;

% sampled velocity components
% each row represents a time instant as dictated by t_smpl
% each column represents a spatial location as dictated by y_smpl
w_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(2));
u_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(3));
v_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(4));


%% instantaneous velocity plots
figure(1)
title('$\textbf{Instantaneous Velocity}$','Interpreter','latex','FontSize',14)
hold on
plot(x_smpl,u_smpl(1,:),':r',LineWidth=2)
plot(x_smpl,v_smpl(1,:),'--b',LineWidth=2)
plot(x_smpl,w_smpl(1,:),'-k',LineWidth=2)
xlabel('$x/\delta$','Interpreter','latex','FontSize',14)
xlim([0,1])
ylim([-0.2,1.4])
legend('$u/u_b$','$v/u_b$','$w/u_b$',...
    'Interpreter','latex','FontSize',14,'Location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',14)
grid on
hold off
%% |Mean Streamwise Velocity Component|
% |Because the normal is mostly regarding the u component, it can be seen that 
% u component is prominent compared to the two other components.| 

%u-velocity sampling
u_mean = mean(u_smpl, 1);
u_std = std(u_smpl);

%v-velocity sampling
v_mean = mean(v_smpl, 1);
v_std = std(v_smpl);


%w-velocity sampling
w_mean = mean(w_smpl, 1);
w_std = std(w_smpl);

%Plotting the Mean and STD
figure(3)
plot(x_smpl,u_mean,':r',LineWidth=2)
hold on
plot(x_smpl,v_mean,'--b',LineWidth=2)
plot(x_smpl,w_mean,'-k',LineWidth=2)
title('$\textbf{Mean Velocity at Every Sampling Point}$','Interpreter','latex','FontSize',14)
xlabel('$y$','Interpreter','latex','FontSize',14)
xlim([0,1])
ylim([-0.2,1.4])
legend('$u/u_b$','$v/u_b$','$w/u_b$',...
    'Interpreter','latex','FontSize',14,'Location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',14)
grid on
hold off

%STD
figure(4)
plot(x_smpl,u_std,':r',LineWidth=2)
hold on
plot(x_smpl,v_std,'--b',LineWidth=2)
plot(x_smpl,w_std,'-k',LineWidth=2)
title('$\textbf{Standard Deviation at Every Sampling Point}$','Interpreter','latex','FontSize',14)
xlabel('$t$','Interpreter','latex','FontSize',14)
legend('$u/u_b$','$v/u_b$','$w/u_b$',...
    'Interpreter','latex','FontSize',14,'Location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',14)
grid on
hold off

%% |Wall Normal Streamwise Velocity Gradient|
% |At this case, only w/ub is the only velocity component that has values in 
% wall-normal direction. The other two has near zero value and thus can be ignored.|

w_grad = gradient(w_mean);
figure(5)
figure(5)
plot(x_smpl, w_grad, '-k', 'LineWidth', 2)
title('$\textbf{Mean Streamwise Velocity Gradient}$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 14)
legend('$w/u_b$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'east')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)

% Get the handle to the title object
title_handle = title('$\textbf{Mean Streamwise Velocity Gradient}$', 'Interpreter', 'latex', 'FontSize', 14);
set(title_handle, 'Position', get(title_handle, 'Position') + [0 0.01 0]); % Adjust the vertical position
grid on
hold off
%% |Analytical|
% |To capture the nearly identical profile to the numerical result, von Kármán 
% log velocity profile is used instead of the usual elliptical profile.| 

%von Kármán Log Velocity Profile
N = 5; %Number of Points
c_f = 0.026/(Re_b)^(1/7);
tau_w = 0.5*c_f*rho*u_b^2;
u_t = sqrt(tau_w/rho);
delta_v = nu/u_t;
kappa = 0.41;
y0 = kappa*delta_v;

syms y
y = linspace(y0, 1, N); % Wall-normal distance array

% Calculate the mean streamwise velocity using the von Kármán logarithmic law
w_a = (u_t / kappa) * log(y / y0);

% Plotting the von Kármán logarithmic law profile
figure(6);
plot(y, w_a, '-o', 'LineWidth', 2);
title('von Kármán Logarithmic Law of the Wall');
xlabel('$x/\delta$','Interpreter','latex','FontSize',14);
grid on
hold off
% |FFT Transformation at x = 0.06 location|
% |Because of the value of the numerical solution array has no exact x = 0.06 
% value, the index of which the value is the closest needs to be determined prior 
% to the calculation| 

%Finding the position of x = 0.06 in the x_smpl array
specified_value = 0.06;

% Find the index of the element in the array closest to the specified value
[~, index] = min(abs(x_smpl - specified_value));

% Display the index and the value in the array that is closest to the specified value
closest_value = x_smpl(index);
disp(['Index of closest element: ', num2str(index)]);
disp(['Value in the array closest to ', num2str(specified_value), ': ', num2str(closest_value)]);

%Plot the noisy signal with time
u_int_area = u_smpl(:, index);
v_int_area = v_smpl(:, index);
w_int_area = w_smpl(:, index);

figure(7)
title('$\textbf{Wall Normal Velocity Signal}$','Interpreter','latex','FontSize',14)
plot(t_smpl,u_int_area,':r',LineWidth=2)
hold on
plot(t_smpl,v_int_area,'--b',LineWidth=2)
plot(t_smpl,w_int_area,'-k',LineWidth=2)
xlabel('$x/\delta$','Interpreter','latex','FontSize',14)
legend('$u/u_b$','$v/u_b$','$w/u_b$',...
    'Interpreter','latex','FontSize',14,'Location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',14)
grid on
hold off
%% 
% |Then the FFT transformation can be made.|


%% 
%