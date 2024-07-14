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
GL = 20;


%CONSTRUCTORES
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);
Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
Rayleigh = classRayleigh(Propagador, sigma, seed);
perturbador = classPerturbator(thermoOpticCoefficient,thermalExpansionCoefficient, poissonRatio, p12, p11, Propagador);
PhiOTDR = classPhiOTDR(10);

Vp = Propagador.Vp;
%%
PhiOTDR.setReference(Rayleigh);

[ReferenceX,ReferenceY] = PhiOTDR.getReference();

T = 3;
e = 6*10^-10;

InitDistance = 23*10^3/segmentLength;
EndDistance = 27*10^3/segmentLength;
perturbador.temperatureChange(T,InitDistance,EndDistance);

InitDistance = 38*10^3/segmentLength;
EndDistance = 39*10^3/segmentLength;
perturbador.strainChange(e,InitDistance,EndDistance);

[changedX,changedY] = Rayleigh.generateRayleighTrace();

[raw,filtered]=PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,GL*10);

fgrMng.newFigure();
plot(Vp*ReferenceX/2000,unwrap(angle(ReferenceY)));
hold on
plot(Vp*ReferenceX/2000,unwrap(angle(changedY)));
hold off
grid on
xlabel("Distancia [Km]")
ylabel("Fase [rad]")
legend("Traza original","Traza perturbada")

fgrMng.newFigure();
plot(Vp*ReferenceX/2000,unwrap(angle(ReferenceY)));
hold on
plot(Vp*ReferenceX/2000,unwrap(angle(changedY)));
hold off
grid on
xlabel("Distancia [Km]")
ylabel("Fase [rad]")
legend("Traza original","Traza perturbada")

fgrMng.newFigure();
plot(Vp*ReferenceX/2000,unwrap(angle(ReferenceY)));
hold on
plot(Vp*ReferenceX/2000,unwrap(angle(changedY)));
hold off
grid on
xlabel("Distancia [Km]")
ylabel("Fase [rad]")
legend("Traza original","Traza perturbada")

fgrMng.newFigure();
plot(Vp*ReferenceX/2000,unwrap(angle(changedY)) - unwrap(angle(ReferenceY)));
grid on
xlabel("Distancia [Km]")
ylabel("Fase [rad]")

fgrMng.newFigure();
plot((1:length(raw))*segmentLength/1000,raw);
grid on
xlabel("Distancia [Km]")
ylabel("Indice de refracción")


fgrMng.newFigure();
plot((1:length(filtered))*segmentLength/1000,filtered);
grid on
xlabel("Distancia [Km]")
ylabel("Indice de refracción")

fgrMng.newFigure();
plot((1:length(filtered))*segmentLength/1000,filtered/(thermoOpticCoefficient+thermalExpansionCoefficient*Propagador.n));
grid on
xlabel("Distancia [Km]")
ylabel("Variación terminca [K]")

fgrMng.newFigure();
epsilon = -(1/2)*(Propagador.n^2)*((1-poissonRatio)*p12 - poissonRatio*p11);
plot((1:length(filtered))*segmentLength/1000,10^9*filtered/(Propagador.n*(epsilon + 1)));
grid on
xlabel("Distancia [Km]")
ylabel("Variación tensorial [n"+[char(949)]+"]")

