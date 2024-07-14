%EFECTO DEL RUIDO
clc
clear
fgrMng = figureManager();

%%

%DATOS DE LA FIBRA
n1 = 1.5;
n2 = 1.46;
coreRadio = 1.7*10^(-6);
Length = 300*10^3;
segmentLength = 1;

%DATOS DEL EMISOR
Wavelength = 1550 *10^-9;
Field = 2*exp(1i*rand*2*pi);
PulseTime = 100*10^-9;


%DATOS DE PROPAGACIÓN
alfaFactor = 1;
betaFactor = 1;

%DATOS DE RAYLEIGH
seed = 3;
sigma = 0.01;

%DATOS DE PERTURBADORES
thermoOpticCoefficient = 3.12*10^(-12);
thermalExpansionCoefficient = 1.72*10^(-12);

poissonRatio=0.17;
p11=0.121;
p12=0.27;

%DATOS GAUGE LENGHT
GL = 30;

%CONSTRUCTORES
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);
Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
Rayleigh = classRayleigh(Propagador, sigma, seed);
perturbador = classPerturbator(thermoOpticCoefficient,thermalExpansionCoefficient, poissonRatio, p12, p11, Propagador);
PhiOTDR = classPhiOTDR(GL);

noiseMedia = 0;
filterWidth = 3*GL;

%%
%Potencia promedio del ruido 10^-6
noiseDesviacion = 10^-6;

%operaciones previas
PhiOTDR.setNoise(noiseMedia,noiseDesviacion);

PhiOTDR.setReference(Rayleigh);
perturbador.temperatureChange(10,100,300)
perturbador.temperatureChange(10,10*10^3,10*10^3+300);
perturbador.temperatureChange(10,60*10^3,60*10^3+300);

perturbador.temperatureChange(10,80*10^3,80*10^3+300);

[raw,filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,filterWidth);
raw = raw/(thermoOpticCoefficient+thermalExpansionCoefficient*Propagador.n);
%filtered = filtered/-(thermoOpticCoefficient+thermalExpansionCoefficient*Propagador.n);


fgrMng.newFigure();
plot((1:length(raw))*segmentLength,raw);
xlabel("Distancia [m]");
ylabel("Variación de temperatura [K]");
grid on;
ylim([-2,15])


fgrMng.newFigure();
plot((1:length(raw))*segmentLength/1000,raw);
xlabel("Distancia [Km]");
ylabel("Variación de temperatura [K]");
grid on;
ylim([-2,15])

fgrMng.newFigure();
plot((1:length(raw))*segmentLength/1000,raw);
xlabel("Distancia [Km]");
ylabel("Variación de temperatura [K]");
grid on;
ylim([-2,15])


fgrMng.newFigure();
plot((1:length(raw))*segmentLength/1000,raw);
xlabel("Distancia [Km]");
ylabel("Variación de temperatura [K]");
grid on;
ylim([-2,15])









