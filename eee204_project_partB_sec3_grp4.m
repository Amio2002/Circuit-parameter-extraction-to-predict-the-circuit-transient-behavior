% Group-4(section-03)
% Part-B
clc;
clear;
data_all = xlsread('PartB_data.xlsx');
ti = data_all(1,:);       
V_ai = data_all(2,:);
V_A(:,1) = V_ai;
Ti(:,1)=ti;

% Constants
R = 3.19073885764119   ;              % Resistance(ohm)
L = 0.0978768658396658;               % Inductance(Henry)
delta_t = ti(2) - ti(1);              % Time step(second)
n = length(ti);                       % Number of time steps
Kt = 0.1;                             % Torque constant
Kb = 0.1;                             % Back EMF constant
bl = 0.15;                            % Load damping coefficient(Ns/rad)
bm = 0.05;                            % Motor damping coefficient(Ns/rad)
bml = bm + bl;                        % Total damping
j = 0.15;                             % Moment of inertia(kgm2)

% aw(i-1)+bw(i)+cw(i+1)=Va
a= ((L*j)/(Kt*delta_t^2)) - (((L*bml)+(R*j))/(2*Kt*delta_t));
b= (((-2*L*j)/(Kt*delta_t^2)) + ((R*bml)/Kt) + Kb);
c=((L*j)/(Kt*delta_t^2)) + (((L*bml)+(R*j))/(2*Kt*delta_t));


% Create matrix A
A = zeros(n, n);
for k = 1:n
    A(k, k) = b;
end
for k = 2:n-1
    A(k, k-1) = a;
end
for k = 1:n-1
    A(k, k+1) = c;
end

% Last condition
A(n, 1) = 1; 
A(n, n) = 0;
V_A(end) =0 ;

% Angular speed omega
omega =inv(A) * V_A;
omega_2=omega(1:n-1);

% Calculated Torque (τ)
for k = 1:n-1
    tau(k) = (j * ((omega(k+1) - omega(k))) / ( delta_t)); 
end
   ti_2=Ti(1:n-1) ; 
   tau_2(:,1)=tau;

% Plot Input Voltage vs Time
figure(1)
plot(Ti, V_A, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Input Voltage V_a (V)');
title('Input Voltage vs. Time');
grid on;

% Angular Speed vs Time
figure(2)
plot(Ti, omega, 'r-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angular Speed \omega (rad/s)');
title('Angular Speed vs. Time');
grid on;

% Torque vs Time
figure(3)
plot(ti_2, tau_2, 'g-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Torque  (N·m)');
title('Torque vs. Time');
grid on;

% Current i_a
i_a = (tau_2+(bml*omega_2))/Kt;
figure(4)
plot(ti_2, i_a, 'm-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Current i_a (A)');
title('Current vs. Time');
grid on;

% Input Power
P_in = V_A(1:n-1) .* i_a; 

% Output Power
P_out = omega_2.*tau_2 ;  

%  Input Power vs Time
figure(5)
plot(ti_2, P_in, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Power (W)');
title('Input Power  vs. Time');
grid on;

% Output Power vs Time 
figure(6)
plot(ti_2, P_out, 'r-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Power (W)');
title('Output Power  vs. Time');
grid on;

% Calculate Efficiency
eta = (P_out ./P_in);  

% Efficiency vs Time
figure(7)
plot(ti_2, eta, 'k-');
xlabel('Time (s)');
ylabel('Efficiency');
title('Efficiency vs. Time');
grid on;