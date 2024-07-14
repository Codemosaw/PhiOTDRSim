%EFECTO DEL RUIDO
clc
clear
fgrMng = figureManager();

%%

%DATOS DE LA FIBRA
n1 = 1.5;
n2 = 1.46;
coreRadio = 1.7*10^(-6);
Length = 1*10^3;
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

%CONSTRUCTORES
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);
Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
Rayleigh = classRayleigh(Propagador, sigma, seed);
perturbador = classPerturbator(thermoOpticCoefficient,thermalExpansionCoefficient, poissonRatio, p12, p11, Propagador);
noise = 10^-5;

fml = Propagador.Vp/(2*fiber_1.fiberLength);
tml = 1/fml;

AN = 2*GL;

Vp = Propagador.Vp;

%%

[t,Eref] = Rayleigh.generateRayleighTrace();
Eref = Eref+noise*randn(1,length(Eref));
Eref = Eref+1i*noise*randn(1,length(Eref));

perturbador.strainChange(3*10^-11,300,500);
perturbador.strainChange(10*10^-11,100,150);
perturbador.strainChange(7*10^-11,600,700);


[t,Echa] = Rayleigh.generateRayleighTrace();
Echa = Echa+noise*randn(1,length(Echa));
Echa = Echa+1i*noise*randn(1,length(Echa));

fgrMng.newFigure();

plot(Vp*t/2,abs(Echa).^2-abs(Eref).^2);
grid on
xlabel("Distancia [m]");
ylabel("Intensidad [a.u.]");


