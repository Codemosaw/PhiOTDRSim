clc
clear

X1 = -10:0.01:10;
X1 = X1 + 1550;
Y1 = custom_normalpdf(X1,1550,2);
X2 = -10:1:10;
X2 = X2 + 1550;
Y2 = custom_normalpdf(X2,1550,2);


figure(1)
plot(X1,Y1);
hold on
stem(X2,Y2);
hold off
grid on
xlabel("Longitud de onda nm")
ylabel("V/m")
xlim
legend("Perfil real","Perfil aproximado")
