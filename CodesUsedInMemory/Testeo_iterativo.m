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


fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);


%DATOS DEL EMISOR
Wavelength = 1550 *10^-9;
Field = 10;
PulseTime = widthSegment*(1/200)*10^-6;


Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);



%DATOS DE PROPAGACIÓN

alfaFactor = 100;%10^14;
betaFactor = 1;%10^(-3)/8;

Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);


[d,t,y]=Propagador.propagateWaveSaveArray(0,fiber_1.numberOfSegments);


%%
%Probar generar traza
sigma = 0.01;
seed = 3;
Rayleigh = classRayleigh(Propagador,sigma, seed);

%%
%nuevas pruebas
%testint the PhiOTDR functionality
Rayleigh.resetAll();
gaugeLength = widthSegment*2+2;
PhiOTDR = classPhiOTDR(gaugeLength);

disp("Generando referencia PhiOTDR...")
PhiOTDR.setReference(Rayleigh);
disp("Referencia PhiOTDR creada!")

[Nx, Ny] = Rayleigh.generateRayleighTrace();
PhaseN = unwrap(angle(Ny));

ip = 500;
fp = 900;

range = (0:2*pi/20:2*pi);
amplitude = 6.3*10^(-9);
range = amplitude*sin(range);

%range = amplitude * triangularSignal(41);

saveMatrix = ones(length(range),length(Ny));

%%
k=1;

Perturbacion_reconstruida = ones(1,length(range));
figure(1)
for deltaN = range
    Rayleigh.resetAll();
    Rayleigh.propagationModule.changeRefractiveIndex(deltaN,ip,fp);
    %disp("Perturbación generada")

    %disp("Generando traza de fases diferenciales...")
    dif = PhiOTDR.getDifferences(Rayleigh);
    %disp("Traza de fases diferenciales generada!")

    [Mx, My] = Rayleigh.generateRayleighTrace();

    PhaseM = unwrap(angle(My));

    subplot(2,2,1);
    plot(PhaseM-PhaseN)
    title("Resta de fases")
    
    ylim([-2.5,2.5])
    %xlim([490,740])
    grid on

    subplot(2,2,2);
    plot(dif/(2*fiber_1.segmentLength*gaugeLength))
    title("Analisis Phi OTDR con gaguge length " + gaugeLength)
    %ylim([-0.1,0.1])
    %xlim([480,740])
    grid on

    subp = subplot(2,2,3);
    plot(PhaseN);
    hold on
    plot(PhaseM);
    title("Fases")
    legend("orginal","modificada")
    hold off
    xlim([ip-100,fp+100])

    grid on
    
    subplot(2,2,4)
    plot(range)
    hold on
    plot(k,range(k),"*r")
    title("Función de cambio")
    hold off
    
    sgtitle("Delta N:" + deltaN + " Delta B:" + deltaN*Propagador.k0 );
    
    %función provisional para encontrar la pendiente
    Perturbacion_reconstruida(k) = sum(dif(500:850))/(-length(dif(500:850))*2*fiber_1.segmentLength*gaugeLength*Propagador.k0); 
    
    k = k + 1;
    
    pause(2)
end

figure(2)
plot(range)
hold on
plot(Perturbacion_reconstruida)
hold off
grid on
legend("Función de perturbación original", "función de perturbación reconstruida")

%%
k=1;

Perturbacion_reconstruida = ones(1,length(range));

for deltaN = range
    Rayleigh.resetAll();
    Rayleigh.propagationModule.changeRefractiveIndex(deltaN,ip,fp);
    %disp("Perturbación generada")

    %disp("Generando traza de fases diferenciales...")
    dif = PhiOTDR.getDifferences(Rayleigh);
    %disp("Traza de fases diferenciales generada!")

    [Mx, My] = Rayleigh.generateRayleighTrace();

    PhaseM = unwrap(angle(My));

    Perturbacion_reconstruida(k) = sum(dif(500:850))/(-length(dif(500:850))*2*fiber_1.segmentLength*gaugeLength*Propagador.k0); 
    
    k = k + 1;
 
end

figure(3)
plot(range)
hold on
plot(Perturbacion_reconstruida)
hold off
grid on
legend("Función de perturbación original", "función de perturbación reconstruida")

