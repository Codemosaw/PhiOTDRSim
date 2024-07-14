clc
clear


%1

%DATOS DE LA FIBRA
n1 = 1.5;
n2 = 1.46;
coreRadio = 1.7*10^(-6);
Length = 500*10^3;
segmentLength = 1;

%DATOS DEL EMISOR
widthSegment = 10;
Wavelength = 1550 *10^-9;
Field = 2;
PulseTime = widthSegment*(1/200)*10^-6;
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);

%DATOS DE PROPAGACIÓN

alfaFactor = 1;%10^14;
betaFactor = 1;%(1/4)*10^(-6);

%DATOS DE RAYLEIGH
seed = 3;
sigma = 0.01;


%figure manager
fgrMng = figureManager();

Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
Rayleigh = classRayleigh(Propagador, sigma, seed);


reflectores = Rayleigh.getReflectorArrays();

fgrMng.newFigure();
histogram(abs(reflectores));
xlabel("amplitud")
title("Reflectores con tamaño "+segmentLength+" m")

fgrMng.newFigure();
histogram(angle(reflectores),'BinWidth',0.1*pi/3);
xlabel("fase")
title("Reflectores con tamaño "+segmentLength+" m")
%%
%500


%DATOS DE LA FIBRA
n1 = 1.5;
n2 = 1.46;
coreRadio = 1.7*10^(-6);
Length = 500*10^3;
segmentLength = 10;

%DATOS DEL EMISOR
widthSegment = 10;
Wavelength = 1550 *10^-9;
Field = 2;
PulseTime = widthSegment*(1/200)*10^-6;
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);

%DATOS DE PROPAGACIÓN

alfaFactor = 1;%10^14;
betaFactor = 1;%(1/4)*10^(-6);

%DATOS DE RAYLEIGH
seed = 3;
sigma = 0.01;


Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
Rayleigh = classRayleigh(Propagador, sigma, seed);


reflectores = Rayleigh.getReflectorArrays();

fgrMng.newFigure();
histogram(abs(reflectores));
xlabel("amplitud")
title("Reflectores con tamaño "+segmentLength+" m")

fgrMng.newFigure();
histogram(angle(reflectores),'BinWidth',0.1*pi/3);
xlabel("fase")
title("Reflectores con tamaño "+segmentLength+" m")
%%
%1Km


%DATOS DE LA FIBRA
n1 = 1.5;
n2 = 1.46;
coreRadio = 1.7*10^(-6);
Length = 500*10^3;
segmentLength = 100;

%DATOS DEL EMISOR
widthSegment = 10;
Wavelength = 1550 *10^-9;
Field = 2;
PulseTime = widthSegment*(1/200)*10^-6;
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);

%DATOS DE PROPAGACIÓN

alfaFactor = 1;%10^14;
betaFactor = 1;%(1/4)*10^(-6);

%DATOS DE RAYLEIGH
seed = 3;
sigma = 0.01;

Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
Rayleigh = classRayleigh(Propagador, sigma, seed);


reflectores = Rayleigh.getReflectorArrays();

fgrMng.newFigure();
histogram(abs(reflectores));
xlabel("amplitud")
title("Reflectores con tamaño "+segmentLength+" m")

fgrMng.newFigure();
histogram(angle(reflectores),'BinWidth',0.1*pi/3);
xlabel("fase")
title("Reflectores con tamaño "+segmentLength+" m")
%%
%5Km


%DATOS DE LA FIBRA
n1 = 1.5;
n2 = 1.46;
coreRadio = 1.7*10^(-6);
Length = 500*10^3;
segmentLength = 10^3;

%DATOS DEL EMISOR
widthSegment = 10;
Wavelength = 1550 *10^-9;
Field = 2;
PulseTime = widthSegment*(1/200)*10^-6;
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);

%DATOS DE PROPAGACIÓN

alfaFactor = 1;%10^14;
betaFactor = 1;%(1/4)*10^(-6);

%DATOS DE RAYLEIGH
seed = 3;
sigma = 0.01;

Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
Rayleigh = classRayleigh(Propagador, sigma, seed);


reflectores = Rayleigh.getReflectorArrays();

fgrMng.newFigure();
histogram(abs(reflectores));
xlabel("amplitud")
title("Reflectores con tamaño "+segmentLength+" m")

fgrMng.newFigure();
histogram(angle(reflectores),'BinWidth',0.1*pi/3);
xlabel("fase")
title("Reflectores con tamaño "+segmentLength+" m")
%%
%10Km


%DATOS DE LA FIBRA
n1 = 1.5;
n2 = 1.46;
coreRadio = 1.7*10^(-6);
Length = 500*10^3;
segmentLength = 10*10^3;

%DATOS DEL EMISOR
widthSegment = 10;
Wavelength = 1550 *10^-9;
Field = 2;
PulseTime = widthSegment*(1/200)*10^-6;
fiber_1 = classFiber(n1, n2, coreRadio, Length, segmentLength);

%DATOS DE PROPAGACIÓN

alfaFactor = 1;%10^14;
betaFactor = 1;%(1/4)*10^(-6);

%DATOS DE RAYLEIGH
seed = 3;
sigma = 0.01;


Emisor_1 = classTransmitter(Wavelength,Field,PulseTime);
Propagador = classPropagation(fiber_1,Emisor_1,alfaFactor,betaFactor);
Rayleigh = classRayleigh(Propagador, sigma, seed);


reflectores = Rayleigh.getReflectorArrays();

fgrMng.newFigure();
histogram(abs(reflectores));
xlabel("amplitud")
title("Reflectores con tamaño "+segmentLength+" m")

fgrMng.newFigure();
histogram(angle(reflectores),'BinWidth',0.2*pi);
xlabel("fase")
title("Reflectores con tamaño "+segmentLength+" m")