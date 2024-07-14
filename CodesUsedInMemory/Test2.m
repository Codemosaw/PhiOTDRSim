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
betaFactor = 1%(1/4)*10^(-6);
disp("Generando propagador...")
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
disp("propagador terminado!")


%Probemos propagar un haz de luz
disp("Generando traza de transimisión...")
[d,t,y_ref]=Propagador.propagateWaveSaveArray(0,fiber_1.numberOfSegments);
disp("Traza de tranmisión terminada!")

thermoOpticCoefficient = 3.12*10^(-10);
thermalExpansionCoefficient = 1.72*10^(-9);

poissonRatio=0.17;
p11=0.121;
p12=0.27;


perturbador = classPerturbator(thermoOpticCoefficient,thermalExpansionCoefficient, poissonRatio, p12, p11, Propagador);

%%
%Variación de temperatura
d_init = 10*10^3;
d_final = 12*10^3;
DeltaNt1 = perturbador.temperatureChange(50,d_init,d_final);

d_init = 100*10^3;
d_final = 120*10^3;
DeltaNt2 = perturbador.temperatureChange(-20,d_init,d_final);

[d,t,y_DeltaT]=Propagador.propagateWaveSaveArray(0,fiber_1.numberOfSegments);

fgrMng.newFigure()

plot(d/1000,(abs(y_ref)-abs(y_DeltaT))*10^12);
grid on
xlabel("Distancia [Km]")
ylabel("Amplitud [pV/m]")


fgrMng.newFigure()

plot(d/1000,unwrap(angle(y_DeltaT))-unwrap(angle(y_ref)));
grid on
xlabel("distancia [Km]")
ylabel("Fase [rad]")


fgrMng.newFigure()

plot(d/1000,unwrap(angle(y_DeltaT))-unwrap(angle(y_ref)));
grid on
xlabel("distancia [Km]")
ylabel("Fase [rad]")

fgrMng.newFigure()

plot(d/1000,unwrap(angle(y_DeltaT))-unwrap(angle(y_ref)));
grid on
xlabel("distancia [Km]")
ylabel("Fase [rad]")

%% Variacion de strain
Propagador.resetPropagationArrays();

d_init = 20*10^3;
d_final = 20.5*10^3;
DeltaNt1 = perturbador.strainChange(0.5*10^-6,d_init,d_final);

d_init = 120*10^3;
d_final = 122*10^3;
DeltaNt2 = perturbador.strainChange(0.9*10^-6,d_init,d_final);

[d,t,y_DeltaT]=Propagador.propagateWaveSaveArray(0,fiber_1.numberOfSegments);

fgrMng.newFigure()

plot(d/1000,((abs(y_ref)-abs(y_DeltaT))*10^12))
grid on
xlabel("distancia [Km]")
ylabel("Amplitud [pV/m]")


fgrMng.newFigure()

plot(d/1000,unwrap(angle(y_DeltaT))-unwrap(angle(y_ref)));
grid on
xlabel("distancia [Km]")
ylabel("Fase [rad]")


fgrMng.newFigure()
plot(d/1000,unwrap(angle(y_DeltaT))-unwrap(angle(y_ref)));
grid on
xlabel("distancia [Km]")
ylabel("Fase [rad]")


fgrMng.newFigure()
plot(d/1000,unwrap(angle(y_DeltaT))-unwrap(angle(y_ref)));
grid on
xlabel("distancia [Km]")
ylabel("Fase [rad]")

%%
Propagador.resetPropagationArrays();
d_init = 90*10^3;
d_final = 120.5*10^3;
DeltaNt1 = perturbador.temperatureChange(-3,d_init,d_final);

d_init = 95*10^3;
d_final = 98*10^3;
DeltaNt2 = perturbador.strainChange(1.1*10^-6,d_init,d_final);

[d,t,y_DeltaT]=Propagador.propagateWaveSaveArray(0,fiber_1.numberOfSegments);

fgrMng.newFigure()

plot(d/1000,(abs(y_ref)-abs(y_DeltaT))*10^12);
grid on
xlabel("distancia [Km]")
ylabel("Amplitud [pV/m]")


fgrMng.newFigure()

plot(d/1000,unwrap(angle(y_DeltaT))-unwrap(angle(y_ref)));
grid on
xlabel("distancia [Km]")
ylabel("Fase [rad]")


fgrMng.newFigure()
plot(d/1000,unwrap(angle(y_DeltaT))-unwrap(angle(y_ref)));
grid on
xlabel("distancia [Km]")
ylabel("Fase [rad]")

fgrMng.newFigure()
plot(d/1000,unwrap(angle(y_DeltaT))-unwrap(angle(y_ref)));
grid on
xlabel("distancia [Km]")
ylabel("Fase [rad]")
