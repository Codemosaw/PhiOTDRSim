clc
clear

%figure manager
fgrMng = figureManager();


%Probar creación de la fibra
%Como referencia valoes usuales para 
%a entre 4 y 6 um
%n1 entre 1.46 a 1.50
%n2 entre 1.44 a 1.46
widthSegment = 10;
%DATOS DE LA FIBRA
n1 = 1.5;
n2 = 1.46;
coreRadio = 1.7*10^(-6);
Length = 200*10^3;
segmentLength = 0.5;

disp("Generando fibra...")
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);
disp("Fibra generada!")

%DATOS DEL EMISOR
Wavelength = 1550 *10^-9;
Field = 2;
PulseTime = widthSegment*(1/200)*10^-6;

disp("Generando transmisor...")
Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
disp("transmisor generado!")


%DATOS DE PROPAGACIÓN

alfaFactor = 1;%10^14;
betaFactor = (1/2)*10^(-6);
disp("Generando propagador...")
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
disp("propagador terminado!")


%Probemos propagar un haz de luz
disp("Generando traza de transimisión...")
[d,t,y]=Propagador.propagateWaveSaveArray(0,fiber_1.numberOfSegments);
disp("Traza de tranmisión terminada!")

fgrMng.newFigure()
subplot(2,1,1)
plot(d/1000,20*log10(abs(y)))
grid on
xlabel("distancia [Km]")
ylabel("Amplitud [dBm]")

subplot(2,1,2)
plot(d,angle(y))
xlim([0,10])
grid on
xlabel("distancia [m]")
ylabel("Fase [rad]")