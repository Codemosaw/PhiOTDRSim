
classdef classPhiOTDR < handle
   properties(Access = private)
       gaugeLength 
       
       referenceTraceX 
       referenceTraceY 
       
       referenceTraceRawX 
       referenceTraceRawY 
       
       referenceDeltaTrace 
    
       noiseU = 0;
       noiseO = 0;
   end
    
   methods
       
       function obj = classPhiOTDR(gaugeLength)
       % brief:
       %    Constructor fot the phiOTDR sensor
       % Sintax:
       %    obj = classPhiOTDR(gaugeLength)
       % inputs:
       %    * gaugeLength - It's the gauge lenght for the calculations
       % outputs:
       %    * obj - instance of the classPhiOTDR class
       % notes:
       %    * To simulate noise in the inputs call the method setNoise
       % 
       % methods:
       %    * setNoise(mu,sigma)
       %    * setGaugeLength(gaugeLength)
       %    * setReference(Rayleigh)
       %    * getReference()
       %    * getReferenceDifferentialPhase()
       %    * getDifferencesOfDiferentialPhase(Rayleigh)
       %    * getDifferencesOfDiferentialPhaseLowPass(Rayleigh,ancho)
       %    * getBetaIndex(rayleigh)
       %    * getDiferentialPhaseTrace(Rayleigh) 
       end
       
       function obj = setNoise(obj,mu,sigma)
       % brief:
       %    Set a noise configurations to be added to the data measurement
       % sintax:
       %   setNoise(mu,sigma)
       % inputs:
       %    * mu - The average noise
       %    * sigma - the standar deviation of the noise
       % returns: NONE
       % notes:
       %    * The noise will be added to both the reference trace and the
       %      trace getted before, made after the call of this method
       %    * It's not necesary to call this method to use the rest of the
       %      methods in the class
       %    * The noise will be simulated as an added gaussian circular
       %      random variable to the fields
            obj.noiseO = sigma;
            obj.noiseU = mu;
           
       end
       
       function obj = setGaugeLength(obj, gaugeLength)
       % brief:
       %    Change gauge lenght 
       % sintax:
       %    setGaugeLength(gaugeLength)
       % inputs:
       %    * gaugeLength - gauge length that will be setted
       % outpus: NONE
       % notes:
       %    * The diferential phase trace will be recalculated whit the
       %      given gauge length
       %    * TODO: it's a good idea to change this function to writeGaugeLength or changeGaugeLength

           obj.gaugeLength = gaugeLength;
           
           for k = gaugeLength+1:length(obj.referenceTraceX)
             obj.referenceDeltaTrace(k-obj.gaugeLength) = obj.referenceTraceY(k) - obj.referenceTraceY(k - obj.gaugeLength); 
          end
       end
       
       
       function obj = setReference(obj,Rayleigh)
          % brief:
          %     Set a reference trace for the DAS/DTS calculations, given
          %     the current conditions of the fiber
          % sintax:
          %     setReference(Rayleigh)
          % inpus:
          %     * Rayleigh - A classRayleigh instance
          % outpus: NONE
          % notes:
          %     * This function MUST be called before any other DAST/DTS 
          %       opperation or changing the envioremental conditions of
          %       the fiber
          %     * The current status of the Rayleigh trace will be used as
          %       reference, the current fiber conditions should be know to
          %       retrieve the actual variations value of the perturbations
          % 
          [obj.referenceTraceX,obj.referenceTraceY] = Rayleigh.generateRayleighTrace(); 
          
          obj.referenceTraceY = obj.referenceTraceY+obj.noiseO*randn(1,length(obj.referenceTraceY)) + obj.noiseU;
          obj.referenceTraceY = obj.referenceTraceY+1i*(obj.noiseO*randn(1,length(obj.referenceTraceY)) + obj.noiseU);
          
          obj.referenceTraceRawX = obj.referenceTraceX;
          obj.referenceTraceRawY = obj.referenceTraceY;
          
          obj.referenceTraceY = unwrap(angle(obj.referenceTraceY));
          
          obj.referenceDeltaTrace = zeros(1,length(obj.referenceTraceX)-obj.gaugeLength);
          
          for k = obj.gaugeLength+1:length(obj.referenceTraceX)
             obj.referenceDeltaTrace(k-obj.gaugeLength) = obj.referenceTraceY(k) - obj.referenceTraceY(k - obj.gaugeLength); 
          end
          
       end
       
       
       function [refX,refY] = getReference(obj)
       % brief:
       %    Get the electric field of the PhiOTDR trace used as reference
       % sintax:
       %    [t,E] = getReference()
       % inputs: NONE
       % outpus:
       %    * t - the time vector in which the reference was measured
       %    * E - the electric field complex vector measured in the PhiOTDR
       %          trace
          refX = obj.referenceTraceRawX;
          refY = obj.referenceTraceRawY;
       end
       
       function diferentialReference = getReferenceDifferentialPhase(obj)
       % brief:
       %    Get the diferential phase trace of the reference PhiOTDR trace
       % sintax:
       %    diferentialReference = getReferenceDifferentialPhase()
       % inputs: NONE
       % outpus:
       %    * diferentialReference - the diferential phase trace of the reference PhiOTDR trace
           diferentialReference = obj.referenceDeltaTrace;
       end
       
       
       function ret = getDifferencesOfDiferentialPhase(obj,Rayleigh)
       % brief:
       %    Get difference of diferential phases traces from a class
       %    classRayleigh instance
       % sintax:
       %    ret = getDifferencesOfDiferentialPhase(Rayleigh)
       % inputs:
       %    * Rayleigh - A class classRayleigh instance
       % outputs:
       %    * ret - The difference of diferential phases traces
          dif = obj.getDiferentialPhaseTrace(Rayleigh);
          ret = dif - obj.referenceDeltaTrace;
       end
       
       
       function [raw,pass] = getDifferencesOfDiferentialPhaseLowPass(obj,Rayleigh,ancho)
       % brief:
       %    Get difference of diferential phases traces from a class
       %    classRayleigh instance filtered whit a low pass
       % sintax:
       %    [raw,pass] = getDifferencesOfDiferentialPhaseLowPass(Rayleigh,ancho)
       % inputs:
       %    * Rayleigh - A class classRayleigh instance
       %    * ancho - The width of the movil average used as low filter
       % outputs:
       %    * raw - The non-filtered difference of diferential phases
       %    * pass - The filtered difference of diferential phases
       % notes:
       %    * As a low pass filter a movil average function is used
           
            raw = obj.getDifferencesOfDiferentialPhase(Rayleigh);
            raw = raw/(-2*Rayleigh.propagationModule.fiber.segmentLength*obj.gaugeLength*Rayleigh.propagationModule.k0); 
            len = length(raw);
            rango = ancho;

            pass = zeros(1,length(rango+1:len));
            %pasa filtros del ancho del gauge lenght
            for k = rango+1:len
               pass(k-rango) = sum(raw(k-rango:k))/(rango+1);
            end
            
       end
       
       
       function [deltaBeta, deltaIndex, pos , dif] = getBetaIndex(obj, rayleigh)
       % brief:
       %    get the beta and index change, the position of the variations
       %    and the difference of diferental phase trace, from a given
       %    instance of a class classRayleigh
       % sintax:
       %    [deltaBeta, deltaIndex, pos , dif] = getBetaIndex(rayleigh)
       % inputs:
       %    * Rayleigh - A class classRayleigh instance
       % outpus:
       %    * deltaBeta - A vector of the variations in the
       %                  propagation parameter beta
       %    * deltaIndex - A vector of the variations in the effective
       %                   refractive index
       %    * pos - a vector of tuples that contains the initial position
       %            of the variation and the final position of the
       %            variations in the first and second elements, measured
       %            in discrete positions
       %    * dif - The difference of diferential phases traces
       % notes:
       %    * This functions tries to obtain the variations given a number
       %      that the variations will be assumed constant
       %    * There isn't any noise in the measurement
       %    * The position tuple is measurement in discrete points, a
       %      interpretations of the results must be made to obtain the range
       %      in meters
       
           dif = obj.getDifferencesOfDiferentialPhase(rayleigh);
           %primero debemos buscar la cantidad de perturbaciones en la
           %fibra, buscamos los puntos en los cuales no hay ceros
           
           %acá se guardaran los puntos iniciales y finales de las
           %perturbaciones
           init_end_matrix = [];
           
           % Inicializar variables
            en_serie = false;
            conteo_series = 0;

            % Iterar sobre el vector
            for i = 2:length(dif)-1
                if (abs(dif(i)) < 10^(-14) ) && (abs(dif(i-1)) < 10^(-14)) && (abs(dif(i+1))< 10^(-14))
                    % Si encontramos un cero y estábamos en una serie, incrementar el conteo
                    if en_serie
                        conteo_series = conteo_series + 1;
                        en_serie = false;
                        endP = i;
                        init_end_matrix = [init_end_matrix ; initP,endP];
                    end
                else
                    % Si encontramos un elemento aleatorio y no estábamos en una serie, marcar el inicio de la serie
                    if ~en_serie
                        en_serie = true;
                        initP = i;
                    end
                end
            end
           
            %Ahora simplemente calculamos el beta y el n en cada segmento
            deltaBeta = zeros(length(init_end_matrix(:,1)),1);
            deltaIndex = zeros(length(init_end_matrix(:,1)),1);
            for filaK = 1:length(init_end_matrix(:,1))
                fila = init_end_matrix(filaK,:);
                initP = fila(1) + floor(obj.gaugeLength/2);
                endP = fila(2) - floor(obj.gaugeLength/2);
                deltaIndex(filaK) = sum(dif(initP:endP))/(-length(dif(initP:endP))*2*rayleigh.propagationModule.fiber.segmentLength*obj.gaugeLength*rayleigh.propagationModule.k0);
                deltaBeta(filaK) = deltaIndex(filaK)*rayleigh.propagationModule.k0;
            end
           pos = init_end_matrix;
       end
       
       function ret = getDiferentialPhaseTrace(obj, Rayleigh) 
         
         [X,Y] = Rayleigh.generateRayleighTrace();
         
         Y = Y + obj.noiseO*randn(1,length(Y))+obj.noiseU;
         Y = Y + 1i*(obj.noiseO*randn(1,length(Y))+obj.noiseU);
         
         Y = unwrap(angle(Y));
         
         ret = zeros(1,length(X)-obj.gaugeLength);
          
          for k = obj.gaugeLength+1:length(X)
             ret(k-obj.gaugeLength) = Y(k) - Y(k - obj.gaugeLength); 
          end
       end
       
   end
   

   %PRIVATES METHODS
   methods(Access = private)
       
   end
   
   
end