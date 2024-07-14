clc
clear

imagen = imread('Modos_de_propagación.png');

figure(1)


V = linspace(0,6.684,1000);

b = -0.0574*V.^4 + 0.2544*V.^3 - 0.2309*V.^2 + 0.1001*V - 0.0036;

plot(V,b,'r', 'LineWidth', 2);
ylim([0,1])
xlim([0,6.684])