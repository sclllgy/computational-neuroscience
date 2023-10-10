%% clean up
close all;
clear;
clc;

%% set parameters
I = 0.5;                  % external stimulus [pA]
C = 1.0;                  % membrane capacitance [μF]
gL =   1.0;  gNa = 10.0;  % membrane conductance [nS]
EL = -70.0;  ENa = 60.0;  % resting or equilibrium potential [mV]

% parameters of steady-state activation (or inactivation) curves
% p_inf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
Vm = -40.0;  Vh = -42.0;
km =  15.0;  kh =  -7.0;

tau_h = 5.0;  % time constant of h_inf [ms]

%% solve transient sodium model
tmin = 0.0;  tmax = 100.0;
interval = [tmin tmax];
X0 = [-50.0, 0.8];
f = @(t, X) transient_sodium(X, I, C, gL, EL, gNa, ENa, Vm, km, Vh, kh, tau_h);
[t1, X1] = ode45(f, interval, X0);

%% caluculate nullclines
xmin = -70.0;  xmax = 50.0;
ymin =   0.0;  ymax =  1.0;
V = linspace(xmin, xmax, 1000)';
[V_null, h_null] = nullcline(V, I, gL, EL, gNa, ENa, Vm, km, Vh, kh);

%% caluculate vector field
X = linspace(xmin, xmax, 30);
Y = linspace(ymin, ymax, 30);
[V_axis, h_axis] = meshgrid(X, Y);
[dVdt, dhdt] = vector_field(V_axis, h_axis, I, C, gL, EL, gNa, ENa, Vm, km, Vh, kh, tau_h);

%% plot
figure(1); hold on;
subplot(2,1,1); hold on;
quiver(V_axis, h_axis, dVdt, dhdt, 4, Marker='.', Alignment='head', ShowArrowHead='off');
plot(V, V_null, 'k', LineWidth=2);
plot(V, h_null, 'k--', LineWidth=2);
plot(X1(:,1), X1(:,2), 'r-', LineWidth=2)
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('membrane voltage, $ V $ [mV]', Interpreter='latex');
ylabel('$$ \rm Na^+ \ inactivation, $$ \it h', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16, YDir='reverse');
grid on;

subplot(2,1,2); hold on;
plot(t1, X1(:,1), LineWidth=2)
xlim([tmin tmax]);
ylim([xmin xmax]);
xlabel('time [ms]', Interpreter='latex');
ylabel('membrane voltage, $ V $ [mV]', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;
