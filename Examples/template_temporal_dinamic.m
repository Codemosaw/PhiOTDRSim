% author: Marcos/Sasmosaw Rojas Mardones
% Este codigo muestra un ejemplo de la utilización
% de la libreria Phase sensitive OTDR en un caso
% estatico, cada linea y seccion del codigo es 
% explicada para mayor claridad de uso

%seccion 0: procesos previos

addpath('lbr'); %añade el directorio a las librerias (clases)

clc  %limpia la terminal
clear %borra todas las variables de memoria


fgrMng = figureManager(); % creacion de un manager de figuras

%%Seccion 1: inicializacion de parametros 

%Datos de la fibra
n1 = 1.5; %indices de refraccion
n2 = 1.46; %indice de refraccion
coreRadio = 1.7*10^(-6); % Radio del nucleo
Length = 2*10^3; %Largo de la fibra en metros
segmentLength = 1; %Largo del segmento en metros

%DATOS DEL EMISOR
Wavelength = 1550 *10^-9; %longitud de onda en metros
Field = 2*exp(1i*rand*2*pi); %fasor del campo electrico
PulseTime = 100*10^-9; %largo del pulso en segundos


%DATOS DE PROPAGACIÓN
alfaFactor = 1;  %escalamiento para la constante de atenuacion
betaFactor = 1;  %escalamiento para la constante de fase

%DATOS DE RAYLEIGH
seed = 3;   %Semilla aleatoria con la que se generaran los distintos reflectores acumulados de Rayleigh
sigma = 0.01; %Desviacion estandar de los reflectores aleatorios acumuladors

%DATOS DE PERTURBADORES
thermoOpticCoefficient = 3.12*10^(-12); %coeficiente termo optico
thermalExpansionCoefficient = 1.72*10^(-12); %coeficiente termo-expansivo

poissonRatio=0.17; %radio de poisson
p11=0.121;  %coeficientes tenso-opticos
p12=0.27;

%DATOS GAUGE LENGHT
GL = 10; %gauge length en NUMERO DE SEGMENTOS

%DATOS RUIDOS
noiseMedia = 0; % valor medio del ruido optico
noiseDesviacion = 0.1*10^-5; %Desviacion del ruido optico

%DATOS FILTRADO
AN = 2*GL;  %Largo del promedio movil utilizado (en caso de usarlo)

%% Seccion 2: Constructores
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

%% Seccion 3: calibracion (en caso de usarla)

PhiOTDR.setReference(Rayleigh); %configura una traza de referencia
DeltaT = 5; %configurar una variacion de temperatura (o strain) que se acerque a los valores que se desean medir
Distancia = 940; %Distancia a la que calibrar 
perturbador.temperatureChange(DeltaT,1,2*10^3); %realizar el cambio de temperatura

gain = 1; %En esta variable se guardara el factor de calibracion final

terminos = 10; %numero de terminos en la media geometrica

%realizacion iterativa de la calibracion
for k = 1:terminos
    [Raw,Filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN); %obtencion de las trazas

    Raw = Raw*FcT; % Conversion a temperaturas
    p = Raw(Distancia); %Temperatura obtenida
    gain = gain*DeltaT/Raw(Distancia); %Correcion del factor de calibracion
end

gain = gain^(1/terminos);

%% Seccion 4: Simulacion

T=0:0.1:20; %Vector temporal

pertTtemp = 10*(1-exp(-linspace(0,10,length(T)))); %Vector de variacion de temperatura


%Las siguientes 3 lineas son prelocación de memoria
[dummyR,dummyF] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN);
saveMatrixFiltered = zeros(length(T),length(dummyF)); %genera matriz para guardar datos filtrados
saveMatrixRaw = zeros(length(T),length(dummyR));    %genera matriz para guardar datos crudos (sin filtrar)


recon = zeros(1,length(T)); %prelocacion de memoria de la reconstruccion de la curva
reconFil = zeros(1,length(T)); %prelocacion de memoria de la reconstruccion de la curva filtrada 
%Variacion dinamica (temporal) de temperatura
for k = 1:length(T)
    100*k/length(T) %imprimir el "porcentaje" de la simulación
    Rayleigh.resetAll(); %deshace cambios previos a la fibra
    
    perturbador.temperatureChange(pertTtemp(k),900,1100); %aplicamos la perturbación termica
    
    [raw,filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN); %recuperamos las trazas
    raw = raw*FcT; %Recuramos la temperatura con el factor de conversion
    filtered = filtered*FcT; 
    saveMatrixFiltered(k,:) = filtered; %guardamos la curva filtrada
    saveMatrixRaw(k,:) = raw; %guardamos la curva corregida con factor de corrección

    recon(k) = raw(Distancia); %reconstruimos la curva (en este caso no se utiliza la funcino filtrada)
    reconFil(k) = filtered(Distancia)
end

%% sección 5: graficar



%Grafico con corrección
fgrMng.newFigure();
plot(T,pertTtemp,'LineWidth',2);
grid on
hold on
plot(T,recon*gain) %Se corrige con el factor de correción
hold off
legend("Perturbación original","Perturbación reconstruida",'LineWidth',2)
xlabel("Tiempo [s]")
ylabel("Variación de temperatura [K]")

%Grafico sin correción
fgrMng.newFigure();
plot(T,pertTtemp,'LineWidth',2);
grid on
hold on
plot(T,recon)
hold off
legend("Perturbación original","Perturbación reconstruida",'LineWidth',2)
xlabel("Tiempo [s]")
ylabel("Variación de temperatura [K]")

%Grafico con filtro
fgrMng.newFigure();
plot(T,pertTtemp,'LineWidth',2);
grid on
hold on
plot(T,reconFil)
hold off
legend("Perturbación original","Perturbación reconstruida",'LineWidth',2)
xlabel("Tiempo [s]")
ylabel("Variación de temperatura [K]")
