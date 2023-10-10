%% clean up
close all;
clear;
clc;

%% set parameters
I =  0.0;  % external stimulus [pA]
a = -0.1;
b =  0.01;
c =  0.02;

%% solve FitzHugh-Nagumo model
tmin = 0.0;  tmax = 300.0;
interval = [tmin tmax];
f = @(t, X) fitzhugh_nagumo(X, I, a, b, c);
X0 = [0.0, 0.2];
[t1, X1] = ode45(f, interval, X0);

%% calculate nullclines
xmin = -0.5;  xmax = 1.2;
ymin = -0.1;  ymax = 0.3;
V = linspace(xmin, xmax, 100);
[V_null, w_null] = nullcline(V, I, a, b, c);

%% calculate vector field
X = linspace(xmin, xmax, 30);
Y = linspace(ymin, ymax, 30);
[V_axis, w_axis] = meshgrid(X, Y);
[dVdt, dwdt] = vector_field(V_axis, w_axis, I, a, b, c);

%% plot
figure(1); hold on;
subplot(2,1,1); hold on;
quiver(V_axis, w_axis, dVdt, dwdt, 4,  Marker='.', Alignment='head', ShowArrowHead='off');
plot(V, V_null, 'k-', LineWidth=2);
plot(V, w_null, 'k--', LineWidth=2);
plot(X1(:,1), X1(:,2), 'r-', LineWidth=2)
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('membrane voltage, $ V $ [mV]', Interpreter='latex');
ylabel('$$ \rm recovery, $$ \it w', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax,'FontSize',16);
grid on;

subplot(2,1,2); hold on;
plot(t1, X1(:,1), '-', LineWidth=2);
xlabel('time [ms]', Interpreter='latex');
ylabel('membrane voltage, $ V $ [mV]', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax,'FontSize',16);
grid on;