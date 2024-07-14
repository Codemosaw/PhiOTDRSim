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

%CONSTRUCTORES
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);
Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
Rayleigh = classRayleigh(Propagador, sigma, seed);
perturbador = classPerturbator(thermoOpticCoefficient,thermalExpansionCoefficient, poissonRatio, p12, p11, Propagador);
PhiOTDR = classPhiOTDR(GL);

noiseMedia = 0;
noiseDesviacion = 1*10^-4;
filterWidth = 3*GL;

PhiOTDR.setNoise(noiseMedia,noiseDesviacion);

fml = Propagador.Vp/(2*fiber_1.fiberLength);
tml = 1/fml;

AN = 2*GL;

%%

PhiOTDR.setReference(Rayleigh);
[dummyR,dummyF] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN); 

%Tiempo
T = 0:1:1000;
T = T*tml;

fpE = 8*10^3;
ipE = 700;
pertEesp = (6*10^-11)*sin(linspace(0,2*pi,50));
pertEtemp = sin(T*2*pi*fpE);

%Prelocación de memoria
saveMatrixFiltered = zeros(length(T),length(dummyF));
saveMatrixRaw = zeros(length(T),length(dummyR));

for k = 1:length(T)
    100*k/1000
    Rayleigh.resetAll();
    %aplicamos la perturbación termica
    %for p = 1:length(pertTesp)
    %    perturbador.temperatureChange(pertTtemp(k)*pertTesp(p),ipT+p,ipT+p);
    %end
    %aplicamos la perturbación espacial
    for p = 1:length(pertEesp)
        perturbador.strainChange(pertEtemp(k)*pertEesp(p),ipE+p,ipE+p);
    end
    [raw,filtered] = PhiOTDR.getDifferencesOfDiferentialPhaseLowPass(Rayleigh,AN);
    
    saveMatrixFiltered(k,:) = filtered;
    saveMatrixRaw(k,:) = raw;
    
end


%%

%Calculo de las FFT

lenPHIOTDRFILTERED =length(saveMatrixFiltered(1,:));
lenPHIOTDRRAW = length(saveMatrixRaw(1,:));
lenT = floor(length(T)/2) + 1;


MatrizFiltered = zeros(lenT,lenPHIOTDRFILTERED);
Matrizraw = zeros(lenT,lenPHIOTDRRAW);


for k = 1:lenPHIOTDRFILTERED

    datos = saveMatrixFiltered(:,k);
    N = length(datos);
    fs = fml; % Cambia según tu frecuencia de muestreo

    % Calcular la FFT de los datos
    fft_resultado = fft(datos);

    % Calcular el vector de frecuencias correspondiente
    frecuenciasFil = (0:N-1) * fs / N;

    % Calcular la magnitud de la FFT (y escalarla)
    magnitud_fft = abs(fft_resultado) / N;

    % Solo graficar la mitad positiva del espectro (simetría)
    mitad_positiva = floor(N/2) + 1;

    MatrizFiltered(:,k) = magnitud_fft(1:mitad_positiva);
    
end



for k = 1:lenPHIOTDRRAW
    datos = saveMatrixRaw(:,k);
    N = length(datos);
    fs = fml; % Cambia según tu frecuencia de muestreo

    % Calcular la FFT de los datos
    fft_resultado = fft(datos);

    % Calcular el vector de frecuencias correspondiente
    frecuenciasRaw = (0:N-1) * fs / N;

    % Calcular la magnitud de la FFT (y escalarla)
    magnitud_fft = abs(fft_resultado) / N;

    % Solo graficar la mitad positiva del espectro (simetría)
    mitad_positiva = floor(N/2) + 1;
    
    Matrizraw(:,k) = magnitud_fft(1:mitad_positiva);
    
end

%%
fgrMng.newFigure();
imagesc([],frecuenciasFil(1:mitad_positiva),2*(MatrizFiltered))
set(gca, 'YDir', 'normal');
%ylim([0,1000])
colorbar

fgrMng.newFigure();
%clims = [0,10];
clim = [0,4.6*10^(-11)];
imagesc([],frecuenciasRaw(1:mitad_positiva),2*(Matrizraw),clim)
set(gca, 'YDir', 'normal');
%ylim([0,1000])
colorbar
%%
fgrMng.newFigure();
frecuenciasFil = frecuenciasFil/1000;

clim = [0,4.6*10^(-11)];
imagesc([],frecuenciasFil(1:mitad_positiva),2*(MatrizFiltered),clim)
set(gca, 'YDir', 'normal');
xlabel("Distancia [m]")
ylabel("Frecuencia [KHz]")
colorbar

fgrMng.newFigure();
clim = [0,4.6*10^(-11)];
imagesc([],frecuenciasFil(1:mitad_positiva),2*(MatrizFiltered),clim)
set(gca, 'YDir', 'normal');
xlabel("puntos")
ylabel("Hz")
colorbar

fgrMng.newFigure();
clim = [0,4.6*10^(-11)];
imagesc([],frecuenciasFil(1:mitad_positiva),2*(MatrizFiltered),clim)
set(gca, 'YDir', 'normal');
xlabel("puntos")
ylabel("Hz")
colorbar

%%
fgrMng.newFigure();
frecuenciasRaw = frecuenciasRaw/1000;

epsilon = -(1/2)*(Propagador.n^2)*((1-poissonRatio)*p12 - poissonRatio*p11);

clim = [0,70];
imagesc([],frecuenciasRaw(1:mitad_positiva),2*(Matrizraw)/(Propagador.n*(epsilon + 1)),clim)
set(gca, 'YDir', 'normal');
xlabel("Distancia [m]")
ylabel("Frecuencia [KHz]")
cb = colorbar(); 
ylabel(cb,'Variacion tensión [\epsilon]','FontSize',10,'Rotation',90)


fgrMng.newFigure();
clim = [0,70];
imagesc([],frecuenciasRaw(1:mitad_positiva),2*(Matrizraw)/(Propagador.n*(epsilon + 1)),clim)
set(gca, 'YDir', 'normal');
xlabel("Distancia [m]")
ylabel("Frecuencia [KHz]")
cb = colorbar(); 
ylabel(cb,'Variacion tensión [\epsilon]','FontSize',10,'Rotation',90)


fgrMng.newFigure();
clim = [0,70];
imagesc([],frecuenciasRaw(1:mitad_positiva),2*10^12*(Matrizraw)/(Propagador.n*(epsilon + 1)),clim)
set(gca, 'YDir', 'normal');
xlabel("Distancia [m]")
ylabel("Frecuencia [KHz]")
cb = colorbar(); 
ylabel(cb,'Variacion tensión [\epsilon]','FontSize',10,'Rotation',90)

