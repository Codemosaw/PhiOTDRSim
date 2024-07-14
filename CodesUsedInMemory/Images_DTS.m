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
Length = 1*10^3;
segmentLength = 1;

disp("Generando fibra...")
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);
disp("Fibra generada!")

%DATOS DEL EMISOR
Wavelength = 1550 *10^-9;
Field = 10;
PulseTime = widthSegment*(1/200)*10^-6;

disp("Generando transmisor...")
Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
disp("transmisor generado!")


%DATOS DE PROPAGACIÓN

alfaFactor = 1;%10^14;
betaFactor = 1;%10^(-3)/8;
disp("Generando propagador...")
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
disp("propagador terminado!")


%Probemos propagar un haz de luz
disp("Generando traza de transimisión...")
[d,t,y]=Propagador.propagateWaveSaveArray(0,fiber_1.numberOfSegments);
disp("Traza de tranmisión terminada!")

sigma = 0.01;
seed = 2;
Rayleigh = classRayleigh(Propagador,sigma, seed);

thermoOpticCoefficient = 10^(-10);
thermalExpansionCoefficient = 10^(-10);

poissonRatio=0.17;
p11=0.121;
p12=0.27;

Vp = Propagador.Vp;

GL = 6;
PhiOTDR = classPhiOTDR(GL);

PhiOTDR.setReference(Rayleigh);

[ReferenceX,ReferenceY] = PhiOTDR.getReference();
fgrMng.newFigure();
plot(angle(ReferenceY));
hold on

Propagador.changeRefractiveIndex(5*10^-9,300,500);
[changedX,changedY] = Rayleigh.generateRayleighTrace();

plot(changedX*Vp/2,angle(changedY));
grid on
xlabel("Distancia [m]")
ylabel("fase [rad]")
legend("Traza original","Traza modificada")
hold off


fgrMng.newFigure();
plot(changedX*Vp/2,unwrap(angle(ReferenceY)));
hold on

plot(changedX*Vp/2,unwrap(angle(changedY)));
grid on
xlabel("Distancia [m]")
ylabel("fase [rad]")
legend("Traza original","Traza modificada")
hold off


referenceDiferential = PhiOTDR.getReferenceDifferentialPhase();
modifiedDiferential = PhiOTDR.getDiferentialPhaseTrace(Rayleigh);

fgrMng.newFigure();
plot((1:length(referenceDiferential))*segmentLength,referenceDiferential);
hold on
plot((1:length(modifiedDiferential))*segmentLength,modifiedDiferential);
grid on
xlabel("Distancia [m]")
ylabel("fase [rad]")
legend("Traza original","Traza modificada")
hold off


fgrMng.newFigure()
plot(modifiedDiferential-referenceDiferential)
grid on
xlabel("Distancia [m]")
ylabel("fase [rad]")



fgrMng.newFigure()
plot((1:length(changedY))*segmentLength,unwrap(angle(changedY))-unwrap(angle(ReferenceY)))
grid on
xlabel("Distancia [m]")
ylabel("fase [rad]")

%%
[raw,filtered]=PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,GL*3);

fgrMng.newFigure();
plot((1:length(raw))*segmentLength,raw)
grid on
xlabel("Distancia [m]")
ylabel("Indice de refracción")

fgrMng.newFigure();
plot((1:length(filtered))*segmentLength,filtered)
grid on
xlabel("Distancia [m]")
ylabel("Indice de refracción")
