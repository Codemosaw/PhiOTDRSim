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
GL = 10;

%CONSTRUCTORES
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);
Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
Rayleigh = classRayleigh(Propagador, sigma, seed);
perturbador = classPerturbator(thermoOpticCoefficient,thermalExpansionCoefficient, poissonRatio, p12, p11, Propagador);
PhiOTDR = classPhiOTDR(10);

Vp = Propagador.Vp;

%%
[refX,refY] = Rayleigh.generateRayleighTrace();
PhiOTDR.setReference(Rayleigh);


ip = 10*10^3/segmentLength;
fp = 30*10^3/segmentLength;

Propagador.changeRefractiveIndex(10^-10,ip,fp);

[chanX,chanY] = Rayleigh.generateRayleighTrace();
[raw,filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,40);

fgrMng.newFigure()

plot(Vp*refX/2000,unwrap(angle(refY)),"b");
hold on
plot(Vp*refX/2000,unwrap(angle(chanY)),"r");
hold on
plot(Vp*refX/2000,unwrap(angle(chanY))-[zeros(1,14282),6.2297*ones(1,length(angle(chanY))-14282)],"r--");

hold off
grid on
xlabel("Distancia [Km]")
ylabel("Fase [rad]")
legend("Traza original", "Traza perturbada con error","Traza perturbada sin error")

fgrMng.newFigure()
plot(Vp*refX/2000,unwrap(angle(chanY))-unwrap(angle(refY)))
xlabel("Distancia [Km]")
ylabel("Fase [rad]")

grid on

fgrMng.newFigure();
plot((1:length(filtered))*segmentLength/1000,filtered)
xlabel("Distancia [Km]")
ylabel("Indice de refracción")
grid on

%%
[refX,refY] = Rayleigh.generateRayleighTrace();
PhiOTDR.setReference(Rayleigh);


ip = 10*10^3/segmentLength;
fp = 30*10^3/segmentLength;

Propagador.changeRefractiveIndex(10^-9,ip,fp);


[raw,filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,40);


fgrMng.newFigure();
plot((1:length(filtered))*segmentLength/1000,filtered)
xlabel("Distancia [Km]")
ylabel("Indice de refracción")
grid on


%%
[refX,refY] = Rayleigh.generateRayleighTrace();
PhiOTDR.setReference(Rayleigh);


ip = 10*10^3/segmentLength;
fp = 30*10^3/segmentLength;

Propagador.changeRefractiveIndex(10^-8,ip,fp);


[raw,filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,40);


fgrMng.newFigure();
plot((1:length(filtered))*segmentLength/1000,filtered)
xlabel("Distancia [Km]")
ylabel("Indice de refracción")
grid on
%%
[refX,refY] = Rayleigh.generateRayleighTrace();
PhiOTDR.setReference(Rayleigh);


ip = 10*10^3/segmentLength;
fp = 30*10^3/segmentLength;

Propagador.changeRefractiveIndex(10^-7,ip,fp);

[raw,filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,40);

fgrMng.newFigure();
plot((1:length(filtered))*segmentLength/1000,filtered)
xlabel("Distancia [Km]")
ylabel("Indice de refracción")
grid on