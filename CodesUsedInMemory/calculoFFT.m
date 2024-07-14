%%Este script es para calcular la fft de perturbaciones especificas
clc

t = 0:10^(-2):2;

perturbacion_1 = 1*10^(-11) * sin(t*4*pi)-1.1*10^(-11);
pert_1_initP = 100;
pert_1_finalP = 300;

perturbacion_2 = 6.3*10^(-10) * triangularSignal(length(t));
pert_2_initP = 500;
pert_2_finalP = 900;

recon_1 = zeros(1,length(t));
recon_2 = zeros(1,length(t));

%vamos a ir añadiendo la perturbación a la fibra
for k = 1:length(t)
    Rayleigh.resetAll();
    
    Rayleigh.propagationModule.changeRefractiveIndex(perturbacion_1(k),pert_1_initP,pert_1_finalP);
    Rayleigh.propagationModule.changeRefractiveIndex(perturbacion_2(k),pert_2_initP,pert_2_finalP);
    
    [deltaBeta, deltaIndex, pos , dif] = PhiOTDR.getBetaIndex(Rayleigh);
    
    recon_1(k) = deltaIndex(1);
    recon_2(k) = deltaIndex(2);
    disp(k+"/"+length(t))
end

figure(1)
subplot(2,1,1)
plot(perturbacion_1)
hold on
plot(recon_1)
hold off
grid on
legend("Función original","función reconstruida")

subplot(2,1,2)
plot(perturbacion_2)
hold on
plot(recon_2)
hold off
grid on
legend("Función original","función reconstruida")


sgtitle("Reconstrucción de la perturbación");

%Calculo de la fft

figure(2)

fs = Rayleigh.fs;
duracion_total = length(perturbacion_2) / fs;

% Calcular la FFT de la señal
fft_1 = fft(recon_1);
fft_2 = fft(recon_2);

% Calcular el vector de frecuencias correspondiente
frecuencias = linspace(0, fs/2, length(fft_1)/2 + 1);

% Calcular el valor absoluto de la parte positiva de la FFT
fft_magnitud_1 = abs(fft_1(1:length(fft_1)/2 + 1));
fft_magnitud_2 = abs(fft_2(1:length(fft_2)/2 + 1));

% Graficar la FFT
subplot(2,1,1)
stem(frecuencias, fft_magnitud_1);
title('FFT de la Señal 1');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
grid on;

subplot(2,1,2)
stem(frecuencias, fft_magnitud_2);
title('FFT de la Señal 2');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
grid on;

