
clc;
clear;
data_all = xlsread ("Sec3_Grp4_partA_fG_data (1)");
f = data_all(1,:);
G = data_all(2,:);
 
% Plotting the Data:
figure(1)
plot( f,G, 'o', 'LineWidth',2,'Color','r'); hold on
xlabel('Frequency(Hz)',fontsize=15);
ylabel('Conductance(1/ohm)',fontsize=15);
 
% Model Fitting:
Y = 1./G.^2;
X = f.^2;
n = length(G);
% Matrix Inversion
A = [n sum(X)
    sum(X) sum(X.^2)];
b = [sum(Y)
    sum(X.*Y)];
a_all = inv(A)*b;
a0 = a_all(1);
a1 = a_all(2);
R = sqrt(a0);
L = sqrt(a1/(4*pi^2));
G_cal=1./(sqrt(R^2+4*pi^2*f.^2*L^2)); 
title('Frequency vs Conductance Data and Fitted Curve', 'fontsize', 12);
% Plot the Fitted Curve:
figure(1)
plot(f,G_cal, 'b', 'LineWidth',2)
legend('Given Data', 'Fitted Curve');
 
% Estimated error
sr = sum((G - G_cal).^2);
% Standard deviation
sd = sqrt(sr./ (n - 2));
 
% Coefficient of determination
y_bar = sum(G)./n;
st = sum((G - y_bar).^2);
r_2 = (st - sr) ./ st;
% Correlation coefficient
r=sqrt(r_2)