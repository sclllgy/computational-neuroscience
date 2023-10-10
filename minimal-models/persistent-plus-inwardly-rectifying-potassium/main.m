%% clean up
close all;
clear;
clc;

%% set parameters
I = 68.0;                 % external stimulus [pA]
C = 1.0;                  % membrane capacitance [μF]
gKir =  20.0;  gK = 2.0;  % membrane conductance [nS]
EK   = -80.0;             % potassium equilibrium potential [mV]

% parameters of steady-state activation curves
% p_inf = 1 ./ (1 + (exp(Vp-V)./kp)), p = h or n
Vh = -80.0;  Vn = -40.0;
kh = -12.0;  kn =   5.0;

tau_n = 5.0;  % time constant of n_inf [ms]

%% solve persistent plus inwardly rectifying potassium model
tmin = 0.0;  tmax = 100.0;
interval = [tmin tmax];
X0 = [-60.0, 0.0];
f = @(t, X) persistent_plus_inwardly_rectifying_potassium(X, I, C, gKir, EK, gK, Vh, kh, Vn, kn, tau_n);
[t1, X1] = ode45(f, interval, X0);

%% caluculate nullcline
xmin = -80.0;  xmax = 30.0;
ymin =  -0.1;  ymax =  1.0;
V = linspace(xmin,xmax,1000)';
[V_null, n_null] = nullcline(V, I, gKir, EK, gK, Vh, kh, Vn, kn);

%% caluculate vector field
X = linspace(xmin, xmax, 30);
Y = linspace(ymin, ymax, 30);
[V_axis, n_axis] = meshgrid(X, Y);
[dVdt, dndt] = vector_field(V_axis, n_axis, I, C, gKir, EK, gK, Vh, kh, Vn, kn, tau_n);

%% plot
figure(1); hold on;
subplot(2,1,1); hold on;
quiver(V_axis, n_axis, dVdt, dndt, 4,  Marker='.', Alignment='head', ShowArrowHead='off');
plot(V, V_null, 'k-', LineWidth=2);
plot(V, n_null, 'k--', LineWidth=2);
plot(X1(:,1), X1(:,2), 'r-', LineWidth=2)
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('membrane voltage, $ V $ [mV]', Interpreter='latex');
ylabel('$$ \rm K^+ \ activation, $$ \it n', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;

subplot(2,1,2); hold on;
plot(t1, X1(:,1), '-', LineWidth=2)
xlim([tmin tmax]);
ylim([xmin xmax]);
xlabel('time [ms]', Interpreter='latex');
ylabel('membrane voltage, $ V $ [mV]', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;