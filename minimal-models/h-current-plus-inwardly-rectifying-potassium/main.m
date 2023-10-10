%% clean up
close all;
clear;
clc;

%% set parameters
I = 10.0;                                % external stimulus [pA]
C = 1.0;                                 % membrane capacitance [μF]
gL =  0.44;  gKir =   4.0;  gh =   0.5;  % membrane conductance [nS]
EL = -50.0;  EK   = -80.0;  Eh = -43.0;  % resting or equilibrium potential [mV]

% parameters of steady-state activation curves
% p_inf = 1./ (1 + (exp(Vp-V)./kp)), p = hKir or h
VhKir = -76.0;  Vh = -65.0;
khKir = -11.0;  kh =  -5.5;

% parameters of voltage-sensitive time constant
% tau_h = C_base + C_amp.*exp(-((V_max-V)./sigma).^2)
C_base =  100.0;
C_amp  = 1000.0;
V_max  =  -75.0;
sigma  =   15.0;

%% solve h-current plus inwardly rectifying potassium
tmin = 0.0;  tmax = 2000.0;
interval = [tmin tmax];
X0 = [-54.0, 0.2];
f = @(t, X) h_current_plus_inwardly_rectifying_potassium(X, I, C, gL, EL, gKir, EK, gh, Eh, VhKir, khKir, Vh, kh, C_base, C_amp, V_max, sigma);
[t1, X1] = ode45(f, interval, X0);

%% caluculate nullcline
xmin = -70.0;  xmax = -43.0;
ymin =   0.0;  ymax =   0.4;
V = linspace(xmin, xmax, 1000)';
[V_null, h_null] = nullcline(V, I, gL, EL, gKir, EK, gh, Eh, VhKir, khKir, Vh, kh);

%% caluculate vector field
X = linspace(xmin, xmax, 30);
Y = linspace(ymin, ymax, 30);
[V_axis, h_axis] = meshgrid(X, Y);
[dVdt, dhdt] = vector_field(V_axis, h_axis, I, C, gL, EL, gKir, EK, gh, Eh, VhKir, khKir, Vh, kh, C_base, C_amp, V_max, sigma);

%% plot
figure(1); hold on;
subplot(2,1,1); hold on;
quiver(V_axis, h_axis, dVdt, dhdt, 4,  Marker='.', Alignment='head', ShowArrowHead='off');
plot(V, V_null, 'k-', LineWidth=2);
plot(V, h_null, 'k--', LineWidth=2);
plot(X1(:,1), X1(:,2), 'r-', LineWidth=2)
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('membrane voltage, $ V $ [mV]', Interpreter='latex');
ylabel('$$ \rm inactivation, $$ \it h', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16, YDir='reverse');
grid on;

subplot(2,1,2); hold on;
plot(t1, X1(:,1), '-', LineWidth=2)
xlim([tmin tmax]);
ylim([xmin xmax]);
xlabel('time [ms]', Interpreter='latex');
ylabel('membrane voltage, $ V $ [mV]', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax,'FontSize',16);
grid on;