clc
clear
fgrMng = figureManager();

%%

%DATOS DE LA FIBRA
n1 = 1.5;
n2 = 1.46;
coreRadio = 1.7*10^(-6);
Length = 10*10^3;
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

%DATOS RUIDO
noiseMedia = 0;
noiseDesviacion = 10^-4;

%CONSTRUCTORES
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);
Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
Rayleigh = classRayleigh(Propagador, sigma, seed);
perturbador = classPerturbator(thermoOpticCoefficient,thermalExpansionCoefficient, poissonRatio, p12, p11, Propagador);
PhiOTDR = classPhiOTDR(GL);

%operaciones previas
PhiOTDR.setNoise(noiseMedia,noiseDesviacion);


%%

fml = Propagador.Vp/(2*fiber_1.fiberLength);
tml = 1/fml;

%%
PhiOTDR.setReference(Rayleigh);

%perfil de cambio
pert = 10*sin(linspace(0,2*pi,1000));

figura = fgrMng.newFigure();

Rayleigh.resetAll();

p = 4000;
n = zeros(1,length(pert));


for indP = pert
    n(p-120+1) = perturbador.temperatureChange(pert(p - 4000 + 1),p,p);
    p = p + 1;
end
[raw,filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,3*GL);

subplot(2,1,1)
plot(raw/-(thermoOpticCoefficient+thermalExpansionCoefficient*Propagador.n));
ylim([-10,10])

subplot(2,1,2)
plot(filtered/-(thermoOpticCoefficient+thermalExpansionCoefficient*Propagador.n));
ylim([-10,10])

pause(0.001)





