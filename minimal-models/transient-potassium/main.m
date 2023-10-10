%% clean up
close all;
clear;
clc;

%% set parameters
I = 10.6;                 % external stimulus [pA]
C = 1.0;                  % membrane capacitance [μF]
gL =   0.2;  gA =   5.0;  % membrane conductance [nS]
EL = -60.0;  EK = -80.0;  % resting or equilibrium potential [mV]

% parameters of steady-state activation curves
% p_inf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
Vm = -45.0;  Vh = -66.0;
km =  10.0;  kh = -10.0;

tau_m = 20.0;  % time constant of m_inf [ms]

%% solve transient potassium (A-current) model
tmin = 0.0;  tmax = 500.0;
interval = [tmin tmax];
X0 = [-25.0, 0.9];
f = @(t, X) transient_potassium(X, I, C, gL, EL, gA, EK, Vm, km, Vh, kh, tau_m);
[t1, X1] = ode45(f, interval, X0);

%% caluculate nullcline
xmin = -80.0;  xmax = 0.0;
ymin =   0.0;  ymax = 1.0;
V = linspace(xmin,xmax,1000)';
[V_null, m_null] = nullcline(V, I, gL, EL, gA, EK, Vm, km, Vh, kh);

%% caluculate vector field
X = linspace(xmin, xmax, 30);
Y = linspace(ymin, ymax, 30);
[V_axis, m_axis] = meshgrid(X, Y);
[dVdt, dmdt] = vector_field(V_axis, m_axis, I, C, gL, EL, gA, EK, Vm, km, Vh, kh, tau_m);

%% plot
figure(1); hold on;
subplot(2,1,1); hold on;
quiver(V_axis, m_axis, dVdt, dmdt, 4,  Marker='.', Alignment='head', ShowArrowHead='off');
plot(V, V_null, 'k-', LineWidth=2);
plot(V, m_null, 'k--', LineWidth=2);
plot(X1(:,1), X1(:,2), 'r-', LineWidth=2)
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('membrane voltage, $ V $ [mV]', Interpreter='latex');
ylabel('$$ \rm K^+ \ activation, $$ \it m', Interpreter='latex');
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