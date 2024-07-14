clc
clear

% Definir los parámetros de la dona
R = 5;  % Radio mayor
r = 5;  % Radio menor

% Definir los ángulos theta y phi
theta = linspace(0, 2*pi, 50);  % Ángulos en el plano horizontal
phi = linspace(0, 2*pi, 50);    % Ángulos en el plano vertical

% Crear las mallas de ángulos theta y phi
[Theta, Phi] = meshgrid(theta, phi);

% Calcular las coordenadas cartesianas para la dona
X = (R + r*cos(Theta)) .* cos(Phi);
Y = (R + r*cos(Theta)) .* sin(Phi);
Z = r*sin(Theta);

% Graficar la dona en 3D
figure(1);

surf(X, Y, Z, 'EdgeColor', 'interp', 'FaceColor', 'interp', 'FaceAlpha', 0);
axis equal;
view(3);
grid off
axis off;