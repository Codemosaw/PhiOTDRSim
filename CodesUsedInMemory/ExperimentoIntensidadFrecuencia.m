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
PhiOTDR = classPhiOTDR(GL);

noiseMedia = 0;
noiseDesviacion = 1*10^-5;
filterWidth = 3*GL;

fml = Propagador.Vp/(2*fiber_1.fiberLength);
tml = 1/fml;

AN = 2*GL;

%%

[t,Eref] = Rayleigh.generateRayleighTrace();
Eref = Eref+noiseDesviacion*randn(1,length(Eref));
Eref = Eref+1i*noiseDesviacion*randn(1,length(Eref));

%Tiempo
TMax = 10000;
T = 0:1:TMax;
T = T*tml;

%Perturbación tensorial
fp1 = 16*10^3;
ip1 = 700;
pertEesp1 = (6*10^-9)*sin(linspace(0,2*pi,49));
pertEtemp1 = sin(T*2*pi*fp1);

fp2 = 20*10^3;
ip2 = 800;
pertEesp2 = (6*10^-9)*sin(linspace(0,3*pi,29));
pertEtemp2 = sin(T*2*pi*fp2);

fp3 = 2*10^3;
ip3 = 100;
pertEesp3 = (6*10^-9)*sin(linspace(0,1.4*pi,59));
pertEtemp3 = sin(T*2*pi*fp3);

%Prelocación de memoria
saveMatrixRaw = zeros(length(T),length(Eref));

for k = 1:length(T)
    k*100/TMax
    %aplicamos la perturbación espacial
    for p = 1:length(pertEesp1)
        perturbador.strainChange(pertEtemp1(k)*pertEesp1(p),ip1+p,ip1+p);
    end
    %aplicamos la perturbación espacial
    for p = 1:length(pertEesp2)
        perturbador.strainChange(pertEtemp2(k)*pertEesp2(p),ip2+p,ip2+p);
    end
    %aplicamos la perturbación espacial
    for p = 1:length(pertEesp3)
        perturbador.strainChange(pertEtemp3(k)*pertEesp3(p),ip3+p,ip3+p);
    end
    
    [t,E] = Rayleigh.generateRayleighTrace();
    E = E + noiseDesviacion*randn(1,length(Eref));
    E = E + 1i*noiseDesviacion*randn(1,length(Eref));
    saveMatrixRaw(k,:) = abs(E).^2-abs(Eref).^2;
    %plot(abs(E).^2-abs(Eref).^2)
    %ylim([-0.1,0.1])
    %pause(1);
end





%%

%Calculo de las FFT

lenPHIOTDRRAW = length(saveMatrixRaw(1,:));
lenT = floor(length(T)/2) + 1;

Matrizraw = zeros(lenT,lenPHIOTDRRAW);


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
frecuenciasRaw = frecuenciasRaw/1000;
fgrMng = figureManager();
fgrMng.newFigure();
clims = [0,1]*10^-3;
imagesc([],frecuenciasRaw(1:mitad_positiva),2*(Matrizraw),clims)
%imagesc([],frecuenciasRaw(1:mitad_positiva),2*(Matrizraw))
set(gca, 'YDir', 'normal');
%ylim([0,1000])
colorbar
xlabel("puntos")
ylabel("Hz")

fgrMng.newFigure();
clims = [0,1]*10^-3;
imagesc([],frecuenciasRaw(1:mitad_positiva),2*(Matrizraw),clims)
%imagesc([],frecuenciasRaw(1:mitad_positiva),2*(Matrizraw))
set(gca, 'YDir', 'normal');
%ylim([0,1000])
colorbar
xlabel("puntos")
ylabel("KHz")

fgrMng.newFigure();
clims = [0,0.05]*10^-3;
imagesc([],frecuenciasRaw(1:mitad_positiva),2*(Matrizraw),clims)
%imagesc([],frecuenciasRaw(1:mitad_positiva),2*(Matrizraw))
set(gca, 'YDir', 'normal');
%ylim([0,1000])
colorbar
xlabel("puntos")
ylabel("KHz")

fgrMng.newFigure();
clims = [0,0.01]*10^-3;
imagesc([],frecuenciasRaw(1:mitad_positiva),2*(Matrizraw),clims)
%imagesc([],frecuenciasRaw(1:mitad_positiva),2*(Matrizraw))
set(gca, 'YDir', 'normal');
%ylim([0,1000])
colorbar
xlabel("puntos")
ylabel("KHz")


