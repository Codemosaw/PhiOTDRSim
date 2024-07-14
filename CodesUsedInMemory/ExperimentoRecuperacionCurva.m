%EFECTO DEL RUIDO
clc
clear
fgrMng = figureManager();

%%

%DATOS DE LA FIBRA
n1 = 1.5;
n2 = 1.46;
coreRadio = 1.7*10^(-6);
Length = 2*10^3;
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
GL = 10;

%Configuraciones de ruido
noiseMedia = 0;
noiseDesviacion = 0.1*10^-5;

%configuraciones de filtrado
filterWidth = 3*GL;

%Ancho para promedio movil
AN = 2*GL;

%CONSTRUCTORES
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);
Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
Rayleigh = classRayleigh(Propagador, sigma, seed);
perturbador = classPerturbator(thermoOpticCoefficient,thermalExpansionCoefficient, poissonRatio, p12, p11, Propagador);
PhiOTDR = classPhiOTDR(GL);
PhiOTDR.setNoise(noiseMedia,noiseDesviacion);

%Frecuencia lenta
fml = Propagador.Vp/(2*fiber_1.fiberLength);
tml = 1/fml;

%Factor de conversión termica
FcT=1/(thermoOpticCoefficient+thermalExpansionCoefficient*Propagador.n);
%Factor de conversion tensorial
epsilon = -(1/2)*(Propagador.n^2)*((1-poissonRatio)*p12 - poissonRatio*p11);
FcE = 1/(Propagador.n*(epsilon + 1));

%% calibración
PhiOTDR.setReference(Rayleigh);
DeltaT = 5; 
Distancia = 940;
perturbador.temperatureChange(DeltaT,1,2*10^3);

gain = 1;

terminos = 10;
for k = 1:terminos
    [Raw,Filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN);

    Raw = Raw*FcT;
    p = Raw(Distancia);
    gain = gain*DeltaT/Raw(Distancia);%media geometrica de la ganancia
end

gain = gain^(1/terminos);

%%

%Rayleigh.resetAll();

T=0:0.1:20;

pertTtemp = 10*(1-exp(-linspace(0,10,length(T))));


[dummyR,dummyF] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN); 

%Prelocación de memoria
saveMatrixFiltered = zeros(length(T),length(dummyF));
saveMatrixRaw = zeros(length(T),length(dummyR));

recon = zeros(1,length(T)); 
for k = 1:length(T)
    100*k/length(T)
    Rayleigh.resetAll();
    %aplicamos la perturbación termica
    perturbador.temperatureChange(pertTtemp(k),900,1100);
    
    [raw,filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN);
    raw = raw*FcT;
    saveMatrixFiltered(k,:) = filtered;
    saveMatrixRaw(k,:) = raw;
    
    %recon(k) = sum(raw(900:1090))/(1091-900);
    recon(k) = raw(Distancia);
end

%%

%El promedio deberia calcularse desde el punto 900 al 1100
%T = T*tml;


fgrMng.newFigure();
plot(T,pertTtemp,'LineWidth',2);
grid on
hold on
plot(T,recon*gain)
hold off
legend("Perturbación original","Perturbación reconstruida",'LineWidth',2)
xlabel("Tiempo [s]")
ylabel("Variación de temperatura [K]")


fgrMng.newFigure();
plot(T,pertTtemp,'LineWidth',2);
grid on
hold on
plot(T,recon)
hold off
legend("Perturbación original","Perturbación reconstruida",'LineWidth',2)
xlabel("Tiempo [s]")
ylabel("Variación de temperatura [K]")


