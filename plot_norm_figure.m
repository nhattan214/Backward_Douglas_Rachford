clc
clear all
close all


% Define the range of x and y values
x = -5:0.01:5; % Adjust the range and resolution as needed
y = -5:0.01:5; % Adjust the range and resolution as needed

% Create a meshgrid for X and Y
[X, Y] = meshgrid(x, y);

% Define the L1 norm function: L1(x, y) = |x| + |y|
L1 = abs(X) + abs(Y)-(sqrt(X.^2+Y.^2));
L2=min(abs(X) + abs(Y),8)-(sqrt(X.^2+Y.^2));
% Plot the contour of L1
figure;
subplot(1,2,1)
contour(X, Y, L1, 20); % 20 contour levels (you can change this number)
colorbar; % Show colorbar for contour levels
xlabel('x');
ylabel('y');
title('$\ell_{1-2}$','Interpreter','latex','FontSize',20);
grid on;
contourcbar("off")

subplot(1,2,2)
contour(X, Y, L2, 20); % 20 contour levels (you can change this number)
colorbar; % Show colorbar for contour levels
xlabel('x');
ylabel('y');
title('$\ell_{1,\mathcal{T}}-\ell_{2}$','Interpreter','latex','FontSize',20);
grid on;
contourcbar("off")
axis tight