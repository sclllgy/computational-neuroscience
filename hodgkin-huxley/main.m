%% clean up
close all;
clear;
clc;

%% set parameters
I1 = 10.0;                             % external stimulus of original Hodgkin-Huxley [pA]
I2 = 5.0;                              % external stimulus of reduced Hodgkin-Huxley [pA]
C = 1.0;                               % membrane capacitance [μF]
gL =  0.3;  gNa = 120.0;  gK =  36.0;  % membrane conductance [nS]
EL = 10.6;  ENa = 120.0;  EK = -12.0;  % resting, equilibrium potential [mV]

%% solve original Hodgkin-Huxley model
tmin = 0.0;  tmax = 50.0;
interval = [tmin tmax];
f = @(t, x) original_hh(x, C, I1, gL, EL, gNa, ENa, gK, EK);
X0_1 = [0.0, 0.1, 0.6, 0.3];
[t1, X1] = ode45(f, interval, X0_1);

%% linear approximation by least squares method
p = polyfit(X1(:,4), X1(:,3), 1);

%% solve reduced Hodgkin-Huxley model
g = @(t, x) reduced_hh(x, C, I2, gL, EL, gNa, ENa, gK, EK, p);
X0_2 = [0.0, 0.0, 0.0, 0.3];
[t2, X2] = ode45(g, interval, X0_2);

%% calculate nullclines of reduced Hodgkin-Huxley model
xmin = -11.0;   xmax = 119;
ymin =   0.0;   ymax = 1.0;
V = linspace(xmin, xmax, 200)';
V_null = zeros(length(V), 1);
n_null = zeros(length(V), 1);
for i = 1:length(V)
    [V_null(i), n_null(i)] = nullcline(V(i), I2, gL, EL, gNa, ENa, gK, EK, p);
end

%% calculate vector field of reduced Hodgkin-Huxley model
X = linspace(xmin, xmax, 30);
Y = linspace(ymin, ymax, 30);
[V_axis, n_axis] = meshgrid(X, Y);
[dVdt, dndt] = vector_field(V_axis, n_axis, C, I2, gL, EL, gNa, ENa, gK, EK, p);

%% plot
figure(1); hold on;
subplot(3,1,1); hold on;
plot(t1, X1(:,1), '-', LineWidth=2);
ylabel('membrane voltage, $ V $ [mV]', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;

subplot(3,1,2); hold on;
plot(t1, X1(:,2), 'r-', LineWidth=2);
plot(t1, X1(:,3), 'g-', LineWidth=2);
plot(t1, X1(:,4), 'b-', LineWidth=2);
ylabel('gating variable', Interpreter='latex');
legend('$m(t)$', '$h(t)$', '$n(t)$', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;

subplot(3,1,3); hold on;
plot(t1, X1(:,3)+X1(:,4), '-', LineWidth=2);
ylim([0.0 2.0]);
xlabel('time [ms]', Interpreter='latex');
ylabel('$h(t) + n(t)$', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;


figure(2); hold on;
subplot(1,1,1); hold on;
plot(X1(:,4), X1(:,3), '.', LineWidth=2);
plot(X1(:,4), p(2)+p(1).*X1(:,4), '-', LineWidth=2);
xlabel('$n(t)$', Interpreter='latex');
ylabel('$h(t)$', Interpreter='latex');
text(0.45, 0.5, ['$h=-$',num2str(abs(p(1))),'$n+$',num2str(abs(p(2)))], Interpreter='latex', FontSize=16);
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;
axis equal


figure(3); hold on;
subplot(1,1,1); hold on;
quiver(V_axis, n_axis, dVdt, dndt, 4, Marker='.', Alignment='head', ShowArrowHead='off');
plot(X2(:,1), X2(:,4), '-', LineWidth=2);
plot(V, V_null, 'k-', LineWidth=2);
plot(V, n_null, 'k--', LineWidth=2);
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('membrane voltage, $ V $ [mV]', Interpreter='latex');
ylabel('$$ \rm K^+ \ activation, $$ \it n', Interpreter='latex');
legend('', 'reduced', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;


figure(4); hold on;
subplot(1,1,1); hold on;
plot(t1, X1(:,1), '-', LineWidth=2);
plot(t2, X2(:,1), '--', LineWidth=2);
xlabel('time [ms]', Interpreter='latex');
ylabel('membrane voltage, $ V $ [mV]', Interpreter='latex');
legend('original','reduced', Interpreter='latex')
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;