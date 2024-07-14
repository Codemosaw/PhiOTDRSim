% author: Marcos/Sasmosaw Rojas Mardones
% Este codigo muestra un ejemplo de la utilización
% de la libreria Phase sensitive OTDR 
% Aqui se revisa:
%   * Perturbacion dinamica temporal y espacialmente
%   * Calibracion distribuida de toda la fibra
%   * Generación de graficos FFT 


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
DeltaE = 7*10^(-12); %configurar una variacion de temperatura (o strain) que se acerque a los valores que se desean medir

perturbador.strainChange(DeltaE,1,2*10^3); %realizar el cambio 

[DummyRaw,DummyFiltered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN); %prelocacion de memoria

gain = ones(1,length(DummyRaw)); %Calcular el largo del vector de ganancias, una ganacia por punto en la traza

terminos = 10; %numero de terminos en la media geometrica

for k = 1:terminos
    [Raw,Filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN); %obtencion de las trazas

    Raw = Raw*FcE;  % Conversion a strain

    gain = (gain*DeltaE)./Raw; %Correcion del factor de calibracion
end

gain = gain.^(1/terminos); %media

%% Seccion 4: Simulacion

% Vector temporal 
% notar que la diferencia con template_temporal_dinamica  es que aqui los puntos se calcularan
% a la maxima frecuencia posible por el simulador
TMax = 10000; %Numero de puntos temporales
T = (0:1:TMax - 1)*tml; %Creacion del vector temporal


%Perturbación tensorial (vibracion)
fp = 16*10^3; %Frecuencia de la perturbacion en Hz
ip = 700; %Punto inicial de la perturbacion en metros
ap = 6*10^-9 %magnitud de la perturbacion en unidades de strain
pertEesp =ap*sin(linspace(0,2*pi,49)); %distribucion espacial de la perturbacion
pertEtemp = sin(T*2*pi*fp);  %Distribucion temporal de la perturbacion


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

    %aplicamos la perturbación espacial
    %basicmanete se aplica la perturbacion punto a punto . . .
    for p = 1:length(pertEesp)
        perturbador.strainChange(pertEtemp(k)*pertEesp(p),ip+p,ip+p); %perturbacion en un unico segmento de la fibra
    end
    
    [raw,filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN); %recuperamos las trazas
    raw = raw*FcT; %Recuramos la temperatura con el factor de conversion
    filtered = filtered*FcT; 
    saveMatrixFiltered(k,:) = filtered; %guardamos la curva filtrada
    saveMatrixRaw(k,:) = raw.*gain; %guardamos la curva corregida con factor de corrección
end
%%Seccion 5: FFT

lenPHIOTDRRAW = length(saveMatrixRaw(1,:)); %Calcular el largo de una traza . . . 
lenT = floor(length(T)/2) + 1; %Calcular el largo para la FFT


MatrizFFTraw = zeros(lenT,lenPHIOTDRRAW); %Prelocacion de memoria para la FFT

for k = 1:lenPHIOTDRRAW
    datos = saveMatrixRaw(:,k); %Extraer todos los puntos temporales en una distancia k
    N = length(datos);  
    fs = fml;   %Calcular la frecuencia de muestreo

    % Calcular la FFT de los datos
    fft_resultado = fft(datos); 

    % Calcular el vector de frecuencias correspondiente
    frecuenciasRaw = (0:N-1) * fs / N;

    % Calcular la magnitud de la FFT (y escalarla)
    magnitud_fft = abs(fft_resultado) / N;

    % Solo graficar la mitad positiva del espectro (simetría)
    mitad_positiva = floor(N/2) + 1;
    
    MatrizFFTraw(:,k) = magnitud_fft(1:mitad_positiva); %Guardar la FFT
    
end


%% sección 6: graficar

%Reconstruccion de la curva en un punto cualquiera, por ejemplo el 570
fgrMng.newFigure();
plot(T,pertEtemp*10^12); %ajustar las unicades a pico
grid on
hold on
plot(T,saveMatrixRaw(:,570)*10^12) %Seleccionar los distintos puntos temporales en la distnacia de 570
hold off
legend("Perturbación original","Perturbación reconstruida")
xlabel("Tiempo [s]")
ylabel("Variación de tensión [p\epsilon]")

%%
fgrMng.newFigure();
frecuenciasRaw = frecuenciasRaw;

clim = [0,10]; %Limite de color...
imagesc([],frecuenciasRaw(1:mitad_positiva),2*abs(MatrizFFTraw)) %Graficar las FFT
set(gca, 'YDir', 'normal');
xlabel("Distancia [m]")
ylabel("Frecuencia [Hz]")
cb = colorbar(); 
ylabel(cb,'Variacion tensión [p\epsilon]','FontSize',10,'Rotation',90)