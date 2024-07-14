clc

gamma_t = 3.12*10^(-12);
eta_t = 1.72*10^(-12);

mu_t=0.17;
p11=0.121;
p12=0.27;

n = Propagador.n;

E = 6*10^-11;
T = 10;

epsilon = -(1/2)*(n^2)*((1-mu_t)*p12 - mu_t*p11);

Le = 1*10^3;
Lt = 4*10^3;

%De = Propagador.k0*n*Le*(epsilon + 1)*E
ne = n*(epsilon + 1)*E
%Dt = Propagador.k0*Lt*(gamma_t+eta_t*n)*T
nt = (gamma_t+eta_t*n)*T

%De+Dt
%%
sigmaS = 0.01;
sigmaR = 10^-5;
alfa = 3.54*10^-5;
L_seg = 1;
E = 2;

p = 300;

SNR = ((sigmaS/sigmaR))*E*exp(-2*alfa*L_seg*p)
%%
sigmaS = 0.01;
sigmaR = 10^-5;
alfa = 3.54*10^-5;
z = 80*10^3;

E = 2;

SNR = ((sigmaS/sigmaR))*E*exp(-2*alfa*z)


%%
Propagador.Vp*2.5*10^(-6)


%%

(3*10^-11)/(thermoOpticCoefficient+thermalExpansionCoefficient*Propagador.n)

%%
epsilon = -(1/2)*(Propagador.n^2)*((1-poissonRatio)*p12 - poissonRatio*p11);

(4.546*10^(-11))/(Propagador.n*(epsilon + 1))