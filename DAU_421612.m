%% |Spatio-Temporal Calculation|
%% |Adapted from read_data.m, generated form Kundu et. al, written by Dr. Tamas Jozsa|
% |Start|

clear
close all
clc
tic
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
mu = nu*rho;

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
% each column represents a spatial location as dictated by x_smpl
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

% Create figure for the combined plot
figure(2);

% First subplot: Mean Velocity at Every Sampling Point
subplot(2,1,1);
plot(x_smpl, u_mean, ':r', 'LineWidth', 2);
hold on;
plot(x_smpl, v_mean, '--b', 'LineWidth', 2);
plot(x_smpl, w_mean, '-k', 'LineWidth', 2);
title('$\textbf{Mean Velocity at Every Sampling Point}$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('$x/\delta$', 'Interpreter', 'latex', 'FontSize', 14);
xlim([0, 1]);
ylim([-0.2, 1.4]);
legend('$u/u_b$', '$v/u_b$', '$w/u_b$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'east');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);
grid on;
hold off;

% Second subplot: Standard Deviation at Every Sampling Point
subplot(2,1,2);
plot(x_smpl, u_std, ':r', 'LineWidth', 2);
hold on;
plot(x_smpl, v_std, '--b', 'LineWidth', 2);
plot(x_smpl, w_std, '-k', 'LineWidth', 2);
title('$\textbf{Standard Deviation at Every Sampling Point}$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('$x/\delta$', 'Interpreter', 'latex', 'FontSize', 14);
xlim([0, 1]);
legend('$u/u_b$', '$v/u_b$', '$w/u_b$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'east');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);
grid on;
hold off;

%% |Wall Normal Mean Streamwise Velocity Gradient|
% |At this case,| $w/u_{b}$ |is the only velocity component that has values 
% in wall-normal direction. The other two has near zero value and thus can be 
% ignored (will have zero gradient in all location)|

dw_dx = gradient(w_mean, x_smpl);


%Analytical Solution Based on symbolic operations
c_f = 0.026/(Re_b)^(1/7);
tau_w = 0.5*c_f*rho*u_b^2;
u_t = 0;
v_t = 0;
w_t = sqrt(tau_w/rho);
delta_v = nu/w_t;
kappa = 0.41;
x0 = kappa*delta_v;


du_dx = gradient(u_mean, x_smpl);
x_smpl_uni = linspace(x0, 1, 126);

syms w(x) x
w = w_t / kappa * log(x / x0);
w_126 = double(subs(w, x, x_smpl_uni));

dw_126 = gradient(w_126, x_smpl_uni);

figure;
subplot(1,2,1);
plot(x_smpl_uni, w_126, 'LineWidth', 2)
hold on
plot(x_smpl, w_mean, 'LineWidth', 2)
title('\textbf{Mean Velocity}', 'Interpreter', 'latex', 'FontSize', 14)
legend('Analytical', 'Numerical', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'east');
xlabel('$x/\delta$', 'Interpreter', 'latex', 'FontSize', 14);
hold off

subplot(1,2,2);
plot(x_smpl_uni, dw_126, 'LineWidth', 2)
hold on
plot(x_smpl, dw_dx, 'LineWidth', 2)
title('\textbf{Mean Velocity Gradient}', 'Interpreter', 'latex', 'FontSize', 14)
legend('Analytical', 'Numerical', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'east');
xlabel('$x/\delta$', 'Interpreter', 'latex', 'FontSize', 14);
hold off

error_1 = norm(dw_dx-dw_126,1)
%%
%ERROR, 252 POINTS
n = length(w_mean);

% Elongate x_smpl
int_x_smpl_252 = NaN(1, 2*n);
int_x_smpl_252(1:2:end) = x_smpl_uni;

for i = 2:2:2*n-2
    int_x_smpl_252(i) = (int_x_smpl_252(i-1) + int_x_smpl_252(i+1))/2;
end

int_x_smpl_252(end) = 0;  % Set the last element to zero
int_x_smpl_252(end) = int_x_smpl_252(end-1);  % Set the last element to the second last element

% Elongate w_mean
int_dw_dx_252 = NaN(1, 2*n);
int_dw_dx_252(1:2:end) = dw_dx;

for i = 2:2:2*n-2
    int_dw_dx_252(i) = (int_dw_dx_252(i-1) + int_dw_dx_252(i+1))/2;
end

int_dw_dx_252(end) = 0;  % Set the last element to zero
int_dw_dx_252(end) = int_dw_dx_252(end-1);  % Set the last element to the second last element

w_252 = double(subs(w, x, int_x_smpl_252));
dw_252 = gradient(w_252, int_x_smpl_252);

dw_252(end) = 0;

figure;
plot(int_x_smpl_252, dw_252)
hold on
plot(int_x_smpl_252, int_dw_dx_252)
hold off

error_2 = norm(int_dw_dx_252-dw_252,1)
%%
%ERROR, 504 POINTS

% Elongate x_smpl
int_x_smpl_504 = NaN(1, 4*n);
int_x_smpl_504(1:2:end) = int_x_smpl_252;

for i = 2:2:4*n-2
    int_x_smpl_504(i) = (int_x_smpl_504(i-1) + int_x_smpl_504(i+1))/2;
end

int_x_smpl_504(end) = 0;  % Set the last element to zero
int_x_smpl_504(end) = int_x_smpl_504(end-1);  % Set the last element to the second last element

% Elongate w_mean
int_dw_dx_504 = NaN(1, 4*n);
int_dw_dx_504(1:2:end) = int_dw_dx_252;

for i = 2:2:4*n-2
    int_dw_dx_504(i) = (int_dw_dx_504(i-1) + int_dw_dx_504(i+1))/2;
end

int_dw_dx_504(end) = 0;  % Set the last element to zero
int_dw_dx_504(end) = int_dw_dx_504(end-1);  % Set the last element to the second last element

w_504 = double(subs(w, x, int_x_smpl_504));
dw_504 = gradient(w_504, int_x_smpl_504);

dw_504(502:504) = 0;

figure;
plot(int_x_smpl_504, dw_504)
hold on
plot(int_x_smpl_504, int_dw_dx_504)
hold off

error_3 = norm(int_dw_dx_504-dw_504,1)
%%
%ERROR, 1008 POINTS

% Elongate x_smpl
int_x_smpl_1008 = NaN(1, 8*n);
int_x_smpl_1008(1:2:end) = int_x_smpl_504;

for i = 2:2:8*n-2
    int_x_smpl_1008(i) = (int_x_smpl_1008(i-1) + int_x_smpl_1008(i+1))/2;
end

int_x_smpl_1008(end) = 0;  % Set the last element to zero
int_x_smpl_1008(end) = int_x_smpl_1008(end-1);  % Set the last element to the second last element

% Elongate w_mean
int_dw_dx_1008 = NaN(1, 8*n);
int_dw_dx_1008(1:2:end) = int_dw_dx_504;

for i = 2:2:8*n-2
    int_dw_dx_1008(i) = (int_dw_dx_1008(i-1) + int_dw_dx_1008(i+1))/2;
end

int_dw_dx_1008(end) = 0;  % Set the last element to zero
int_dw_dx_1008(end) = int_dw_dx_1008(end-1);  % Set the last element to the second last element

w_1008 = double(subs(w, x, int_x_smpl_1008));
dw_1008 = gradient(w_1008, int_x_smpl_1008);

dw_1008(1002:1008) = 0;

figure;
plot(int_x_smpl_1008, dw_1008)
hold on
plot(int_x_smpl_1008, int_dw_dx_1008)
hold off

error_4 = norm(int_dw_dx_1008-dw_1008,1)
%%
%ERROR, 2016 POINTS

% Elongate x_smpl
int_x_smpl_2016 = NaN(1, 16*n);
int_x_smpl_2016(1:2:end) = int_x_smpl_1008;

for i = 2:2:16*n-2
    int_x_smpl_2016(i) = (int_x_smpl_2016(i-1) + int_x_smpl_2016(i+1))/2;
end

int_x_smpl_2016(end) = 0;  % Set the last element to zero
int_x_smpl_2016(end) = int_x_smpl_2016(end-1);  % Set the last element to the second last element

% Elongate w_mean
int_dw_dx_2016 = NaN(1, 16*n);
int_dw_dx_2016(1:2:end) = int_dw_dx_1008;

for i = 2:2:16*n-2
    int_dw_dx_2016(i) = (int_dw_dx_2016(i-1) + int_dw_dx_2016(i+1))/2;
end

int_dw_dx_2016(end) = 0;  % Set the last element to zero
int_dw_dx_2016(end) = int_dw_dx_2016(end-1);  % Set the last element to the second last element

w_2016 = double(subs(w, x, int_x_smpl_2016));
dw_2016 = gradient(w_2016, int_x_smpl_2016);

dw_2016(2002:2016) = 0;

figure;
plot(int_x_smpl_2016, dw_2016)
hold on
plot(int_x_smpl_2016, int_dw_dx_2016)
hold off

error_5 = norm(int_dw_dx_2016-dw_2016,1)
%%
error_total = [error_1 error_2 error_3 error_4 error_5]
mesh_size = [126 252 504 1008 2016]

figure;
plot(mesh_size, error_total)

p = polyfit(mesh_size(2:end), error(2:end), 1);  % fit a straight line to the log-log plot
error_Richardson = abs(4*(error(2:end) - error(1:end-1))/ 3);

error_comb = [error_total error_Richardson]
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

figure(6)
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

% Compute FFT 
Fs = 1 / (t_smpl(2) - t_smpl(1));  % Sampling frequency
L = length(u_int_area);  % Length of the signal

%Computing standard FFT
u_fft = fft(u_int_area);
P2_u_fft = abs(u_fft / L);
P1_u_fft = P2_u_fft(1:L/2+1);
f_u_fft = Fs*(0:(L/2))/L;

v_fft = fft(v_int_area);
P2_v_fft = abs(v_fft / L);
P1_v_fft = P2_v_fft(1:L/2+1);
f_v_fft = Fs*(0:(L/2))/L;

w_fft = fft(w_int_area);
P2_w_fft = abs(w_fft / L);
P1_w_fft = P2_w_fft(1:L/2+1);
f_w_fft = Fs*(0:(L/2))/L;

figure(7)
fftlin_1 = subplot (3,3,1);
plot(f_u_fft, P1_u_fft, 'r', 'LineWidth', 2)
title('\textbf{$u/u_b$}','Interpreter','latex','FontSize',14)
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

fftlin_2 = subplot(3,3,2);
plot(f_v_fft, P1_v_fft, 'b', 'LineWidth', 2)
title('\textbf{$v/u_b$}', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

fftlin_3 = subplot(3,3,3);
plot(f_u_fft, P1_u_fft, 'k', 'LineWidth', 2)
title('\textbf{$w/u_b$}', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;



fftsl_1 = subplot (3,3,4);
semilogx(f_u_fft, P1_u_fft, 'r', 'LineWidth', 2)
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

fftsl_2 = subplot(3,3,5);
semilogx(f_v_fft, P1_v_fft, 'b', 'LineWidth', 2)
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

fftsl_3 = subplot(3,3,6);
semilogx(f_w_fft, P1_w_fft, 'k', 'LineWidth', 2)
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;


fftlog_1 = subplot(3,3,7);
loglog(f_u_fft, P1_u_fft, 'r', 'LineWidth', 2)
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

fftlog_2 = subplot(3,3,8);
loglog(f_v_fft, P1_v_fft, 'b', 'LineWidth', 2)
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

fftlog_3 = subplot(3,3,9);
loglog(f_u_fft, P1_u_fft, 'k', 'LineWidth', 2)
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

hold off


% Apply window functions to the signals
window_types = {'rectwin', 'hann', 'hamming', 'blackman'}; % Window types to try

figure;
for i = 1:length(window_types)
    windowed_u = u_int_area .* feval(window_types{i}, length(u_int_area)); % Apply window function to u_int_area
    windowed_v = v_int_area .* feval(window_types{i}, length(v_int_area)); % Apply window function to v_int_area
    windowed_w = w_int_area .* feval(window_types{i}, length(w_int_area)); % Apply window function to w_int_area
    
    % Compute single-sided FFT for windowed signals
    Y_u = fft(windowed_u);
    P2_u = abs(Y_u / L);
    P1_u = P2_u(1:L/2+1);
    f_u = Fs*(0:(L/2))/L;

    Y_v = fft(windowed_v);
    P2_v = abs(Y_v / L);
    P1_v = P2_v(1:L/2+1);
    f_v = Fs*(0:(L/2))/L;

    Y_w = fft(windowed_w);
    P2_w = abs(Y_w / L);
    P1_w = P2_w(1:L/2+1);
    f_w = Fs*(0:(L/2))/L;

    subplot(length(window_types),3,3*i-2);
    loglog(f_u, P1_u, 'r', 'LineWidth', 2);
    title(['$u/u_b$ signal with ', window_types{i}, ' Window'], 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;

    subplot(length(window_types),3,3*i-1);
    loglog(f_v, P1_v, 'b', 'LineWidth', 2);
    title(['$v/u_b$ signal with ', window_types{i}, ' Window'], 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;

    subplot(length(window_types),3,3*i);
    loglog(f_w, P1_w, 'k', 'LineWidth', 2);
    title(['$w/u_b$ signal with ', window_types{i}, ' Window'], 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
end
%%
% Define parameters for filtering
filter_order = 4;  % Butterworth filter order
cutoff_frequency = 0.1;  % Butterworth filter cutoff frequency normalized to Nyquist frequency

% Perform FFT and apply different filters for u, v, and w
u_fft_ma = movmean(u_fft, 5);  % Moving average filtering for u
v_fft_ma = movmean(v_fft, 5);  % Moving average filtering for v
w_fft_ma = movmean(w_fft, 5);  % Moving average filtering for w

u_fft_gaussian = imgaussfilt(abs(u_fft));  % Gaussian filtering for u
v_fft_gaussian = imgaussfilt(abs(v_fft));  % Gaussian filtering for v
w_fft_gaussian = imgaussfilt(abs(w_fft));  % Gaussian filtering for w

[b, a] = butter(filter_order, cutoff_frequency);  % Butterworth filter coefficients
u_fft_butter = filter(b, a, u_fft);  % Butterworth filtering for u
v_fft_butter = filter(b, a, v_fft);  % Butterworth filtering for v
w_fft_butter = filter(b, a, w_fft);  % Butterworth filtering for w

% FFT Filtering Results for Signal u
subplot(3, 3, 1);
loglog(f_u_fft, abs(u_fft_ma(1:L/2+1)), 'r');
title('$u/u_b$ (Moving Average)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Magnitude', 'Interpreter', 'latex', 'FontSize', 14);

subplot(3, 3, 4);
loglog(f_u_fft, abs(u_fft_gaussian(1:L/2+1)), 'r');
title('$u/u_b$ (Gaussian)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Magnitude', 'Interpreter', 'latex', 'FontSize', 14);

subplot(3, 3, 7);
loglog(f_u_fft, abs(u_fft_butter(1:L/2+1)), 'r');
title('$u/u_b$ (Butterworth)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Magnitude', 'Interpreter', 'latex', 'FontSize', 14);

% FFT Filtering Results for Signal v
subplot(3, 3, 2);
loglog(f_v_fft, abs(v_fft_ma(1:L/2+1)), 'b');
title('$v/u_b$ (Moving Average)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Magnitude', 'Interpreter', 'latex', 'FontSize', 14);

subplot(3, 3, 5);
loglog(f_v_fft, abs(v_fft_gaussian(1:L/2+1)), 'b');
title('$v/u_b$ (Gaussian)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Magnitude', 'Interpreter', 'latex', 'FontSize', 14);

subplot(3, 3, 8);
loglog(f_v_fft, abs(v_fft_butter(1:L/2+1)), 'b');
title('$v/u_b$ (Butterworth)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Magnitude', 'Interpreter', 'latex', 'FontSize', 14);

% FFT Filtering Results for Signal w
subplot(3, 3, 3);
loglog(f_w_fft, abs(w_fft_ma(1:L/2+1)), 'k');
title('$w/u_b$ (Moving Average)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Magnitude', 'Interpreter', 'latex', 'FontSize', 14);

subplot(3, 3, 6);
loglog(f_w_fft, abs(w_fft_gaussian(1:L/2+1)), 'k');
title('$w/u_b$ (Gaussian)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Magnitude', 'Interpreter', 'latex', 'FontSize', 14);

subplot(3, 3, 9);
loglog(f_w_fft, abs(w_fft_butter(1:L/2+1)), 'k');
title('$w/u_b$ (Butterworth)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Magnitude', 'Interpreter', 'latex', 'FontSize', 14);


%% |Skewness and Kurtosis|

usmplvect = u_smpl(:);
vsmplvect = v_smpl(:);
wsmplvect = w_smpl(:);

% Normality tests
fprintf('Skewness of u_smpl: %.4f\n', skewness(usmplvect(:)));
fprintf('Kurtosis of u_smpl: %.4f\n', kurtosis(usmplvect(:)));

fprintf('Skewness of v_smpl: %.4f\n', skewness(vsmplvect(:)));
fprintf('Kurtosis of v_smpl: %.4f\n', kurtosis(vsmplvect(:)));

fprintf('Skewness of w_smpl: %.4f\n', skewness(wsmplvect(:)));
fprintf('Kurtosis of w_smpl: %.4f\n', kurtosis(wsmplvect(:)));

% Calculate the histogram
[countsu, edgesu] = histcounts(usmplvect);
[countsv, edgesv] = histcounts(vsmplvect);
[countsw, edgesw] = histcounts(wsmplvect);

% Define the range of x values
x_values_u = edgesu(1:end-1);
x_values_v = edgesv(1:end-1);
x_values_w = edgesw(1:end-1);

% Create a figure with a 3-by-1 grid of subplots
figure;
sgtitle('\textbf{Probability Density Function}','Interpreter', 'latex', 'FontSize', 14)
% Plot the PDF of u_smpl in the first subplot
ax1 = subplot(3,1,1);
plot(x_values_u, countsu);
ylabel('Probability');
title('$u/u_b$', 'Interpreter', 'latex', 'FontSize', 12);

% Plot the PDF of v_smpl in the second subplot
ax2 = subplot(3,1,2);
plot(x_values_v, countsv);
ylabel('Probability');
title('$v/u_b$', 'Interpreter', 'latex', 'FontSize', 12');

% Plot the PDF of w_smpl in the third subplot
ax3 = subplot(3,1,3);
plot(x_values_w, countsw);
ylabel('Probability');
title('$w/u_b$', 'Interpreter', 'latex', 'FontSize', 12');


%%
% Correlation analysis
alpha = 0.05; % Significance level
[r_uv, p_uv] = corrcoef(u_smpl(:), v_smpl(:));
[r_uw, p_uw] = corrcoef(u_smpl(:), w_smpl(:));
[r_vw, p_vw] = corrcoef(v_smpl(:), w_smpl(:));

fprintf('Correlation between u_smpl and v_smpl: %.4f \n', r_uv(1,2));
fprintf('Correlation between u_smpl and w_smpl: %.4f \n', r_uw(1,2));
fprintf('Correlation between v_smpl and w_smpl: %.4f \n', r_vw(1,2));

toc