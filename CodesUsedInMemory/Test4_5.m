clc
clear
fgrMng = figureManager();

%%

%DATOS DE LA FIBRA
n1 = 1.5;
n2 = 1.46;
coreRadio = 1.7*10^(-6);
Length = 50*10^3;
segmentLength = 1;

%DATOS DEL EMISOR
Wavelength = 1550 *10^-9;
Field = 2*exp(1i*rand);
PulseTime = 100*10^-9;
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);

%DATOS DE PROPAGACIÓN

alfaFactor = 0;%10^14;
betaFactor = 1;%(1/4)*10^(-6);

%DATOS DE RAYLEIGH
seed = 3;
sigma = 0.01;


Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
Rayleigh = classRayleigh(Propagador, sigma, seed);

[t,E] = Rayleigh.generateRayleighTrace();

fgrMng.newFigure();
plot(t,abs(E));
grid on
xlabel("tiempo s")
ylabel("V/m")

fgrMng.newFigure();
plot(t,angle(E));
grid on
xlabel("tiempo s")
ylabel("rad")

fgrMng.newFigure();
plot(t,unwrap(angle(E)));
grid on
xlabel("tiempo s")
ylabel("rad")

fgrMng.newFigure();
histogram(abs(E));
grid on
xlabel("Amplitud [V/m]")
ylabel("frecuencia")

fgrMng.newFigure();
histogram(angle(E),'BinWidth',0.1*pi/3);
grid on
xlabel("rad")

fgrMng.newFigure();
plot(abs(E(1:100:end)).^2)
grid on
xlabel("tiempo s")
ylabel("Magnitud")
