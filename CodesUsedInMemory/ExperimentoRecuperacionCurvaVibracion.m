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
noiseDesviacion = (1/4)*10^-5;

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


%Factor de conversion tensorial
epsilon = -(1/2)*(Propagador.n^2)*((1-poissonRatio)*p12 - poissonRatio*p11);
FcE = 1/(Propagador.n*(epsilon + 1));

%% Calibración

PhiOTDR.setReference(Rayleigh);
DeltaE = 10*10^(-12); 
Distancia = 940;
perturbador.strainChange(DeltaE,1,2*10^3);

[DummyRaw,DummyFiltered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN);

gain = ones(1,length(DummyRaw));

terminos = 30;
for k = 1:terminos
    [Raw,Filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN);

    Raw = Raw*FcE;

    gain = (gain*DeltaE)./Raw;
end

gain = gain.^(1/terminos);

%% perturbacion
%PhiOTDR.setReference(Rayleigh);

TMax = 500;

T = 0:1:TMax;
T = T*tml;


fp1 = 500;
ipE = 500;
pertEtemp = 7*10^(-12)*sin(T*2*pi*fp1);
%pertEesp =  sin(linspace(0,2*pi,100));

[dummyR,dummyF] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN); 

%Prelocación de memoria
saveMatrixFiltered = zeros(length(T),length(dummyF));
saveMatrixRaw = zeros(length(T),length(dummyR));


%APLICACION PERTURBACION TEMPORALMENTE
for k = 1:length(T)
    100*k/length(T)
    Rayleigh.resetAll();
    
    
    %APLICACION PERTURBACION ESPACIALMENTE
    
    perturbador.strainChange(pertEtemp(k),ipE,ipE+150);
    
    
    %OBTENCIÓN DE DELTA N
    [raw,filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN);
    
    raw = raw*FcE;
    
    saveMatrixFiltered(k,:) = filtered;
    saveMatrixRaw(k,:) = raw.*gain;
    
end
%% FFT

lenPHIOTDRRAW = length(saveMatrixRaw(1,:));
lenT = floor(length(T)/2) + 1;


MatrizFFTraw = zeros(lenT,lenPHIOTDRRAW);

for k = 1:lenPHIOTDRRAW
    datos = saveMatrixRaw(:,k);
    N = length(datos);
    fs = fml;

    % Calcular la FFT de los datos
    fft_resultado = fft(datos);

    % Calcular el vector de frecuencias correspondiente
    frecuenciasRaw = (0:N-1) * fs / N;

    % Calcular la magnitud de la FFT (y escalarla)
    magnitud_fft = abs(fft_resultado) / N;

    % Solo graficar la mitad positiva del espectro (simetría)
    mitad_positiva = floor(N/2) + 1;
    
    MatrizFFTraw(:,k) = magnitud_fft(1:mitad_positiva);
    
end



%%

%El promedio deberia calcularse desde el punto 900 al 1100
%T = T*tml;

fgrMng.newFigure();
plot(T*1000,pertEtemp*10^12);
grid on
hold on
plot(T*1000,saveMatrixRaw(:,570)*10^12)
hold off
legend("Perturbación original","Perturbación reconstruida")
xlabel("Tiempo [ms]")
ylabel("Variación de tensión [p\epsilon]")

%%
fgrMng.newFigure();
frecuenciasRaw = frecuenciasRaw;

clim = [0,7];
imagesc([],frecuenciasRaw(1:mitad_positiva),2*(MatrizFFTraw)*10^12,clim)
set(gca, 'YDir', 'normal');
xlabel("Distancia [m]")
ylabel("Frecuencia [Hz]")
cb = colorbar(); 
ylabel(cb,'Variacion tensión [p\epsilon]','FontSize',10,'Rotation',90)
