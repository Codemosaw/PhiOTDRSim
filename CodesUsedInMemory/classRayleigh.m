% author: Marcos/Sasmosaw Rojas Mardonesclassdef classRayleigh < handle
   properties
       fs
       
       propagationModule
   end
   
   properties(Access = private)
       
       seed
       sigma
       
       reflectorArray
       reflectorArrayBackup

   end
    
   methods
       
       
       function obj = classRayleigh(propagationModule, sigma, seed)
       % brief: 
       %        The constructor for the Rayleigh class, the
       %        Rayleigh class contain all the data and methods to
       %        generate phase sensitive traces
       % sintax:
       %        classRayleigh(propagationModule, sigma, seed)
       % inputs:
       %    * propagationModule - A classPropagation  instance,
       %                          all the data for propagation will
       %                          be acquired from this instance
       %    * sigma - standar deviaton of the gaussian
       %              related whit how many and how big are the
       %              imperfections in the fiber, a value 0 from
       %              sigma means an ideal fiber.
       %    * seed - Seed for the random generation of reflector
       % outputs:
       %    instance of the classRayleigh class
       %
       % methods:
       %    * resetReflectorArrays()
       %    * resetAll()
       %    * generateRayleighTrace()
       %    * changeReflectors(changePhi, initPoint, finalPoint)
       %    * getReflectorArrays()
       %
       % note:
       %    * Here's (in this class) a more or less good example 
       %      of how the privacy should be manage in a class
               obj.propagationModule = propagationModule;
               obj.seed = seed;
               obj.sigma = sigma;
               rng(seed)
               
               numberOfSegments = obj.propagationModule.fiber.numberOfSegments;
               
               
               obj.reflectorArray = sigma*randn(1,numberOfSegments) + 1i*sigma*randn(1,numberOfSegments);
               obj.reflectorArrayBackup = obj.reflectorArray;
               %la frecuencia de muestreo máximo viene dada por el tiempo
               %que tarda una onda de luz en viajar desde z=0 hasta el
               %final de la fibra, es decir
               obj.fs = propagationModule.Vp/(propagationModule.fiber.segmentLength*2);
       end
           
       
       
       function obj = resetReflectorArrays(obj)
       % brief:
       %    Reset the reflector array
       % sintax:
       %    resetReflectorArrays()
       % inputs: NONE
       % outpus: NONE
           obj.reflectorArray = obj.reflectorArrayBackup;
       end
       
       function obj = resetAll(obj)
       % brief:
       %    Reset all the data modules
       % sintax:
       %    resetAll()
       % inputs: NONE
       % outpus: NONE
            obj.resetReflectorArrays();
            obj.propagationModule.resetPropagationArrays(); 
       end
       
       
       function [t,E] = generateRayleighTrace(obj)
       % brief:
       %        Generate a Phase Sensitive OTDR trace or Rayleigh
       %        trace
       % sintax:
       %        [t,E] = generateRayleighTrace()
       % Description:
       %        This method use all the current instances to get the
       %        data of the fiber, emitter and propagation to produce a
       %        Phase Sensitive OTDR trace whit the vector of times
       %        associetted whit each measurament.
       % inputs:
       %    None
       % outputs:
       %    * t - It's the time vector associetted whit the measurement
       %          of the photodetector
       %    * E - It's the field vector of each measurament made in the
       %          trace
       %
           
           
           numberOfSegments = obj.propagationModule.fiber.numberOfSegments;
           pulseWidth = floor(obj.propagationModule.widthSegments/2);
           [Distance,Time,Field]=obj.propagationModule.propagateWaveSaveArray(0,numberOfSegments);
           Field = Field(1:end-1);
           limit = numberOfSegments - pulseWidth + 1;
           E = zeros(1, limit);
           t = 0:length(E)-1;
           
           for k = 1:limit
               pulseInitPos = k - 1;
               %factor = abs(Field(k))*exp(2*1i*angle(Field(k)));
               factor = Field(k)^2;
               
               suma = 0;
               for finalSegmento = pulseInitPos + 1: pulseInitPos + pulseWidth
                   
                   fase = obj.propagateSingleWave(1,pulseInitPos,finalSegmento);
                   %fase = 2*angle(fase);
                   fase = fase^2;
                   %suma = suma + obj.reflectorArray(finalSegmento)*exp(1i*fase);
                   suma = suma + obj.reflectorArray(finalSegmento)*fase;
               end
               E(k) = factor * suma;
           end
           
           t = t/obj.fs;
            
           
       end
       
       
      
       
       
       
       function obj = changeReflectors(obj,changePhi, initPoint, finalPoint)
       % @brief Change the reflector number to represent changes in
       % refraction index
       % @param in change Factor change, its a complex number where the
       % final reflectore will be given by
       % Rf = Ri * change
       % @param in initPoint Initial point where the change will begin
       % @param in finalPoint Final point where the change will end
           for k = initPoint:finalPoint
              obj.reflectorArray(k) = obj.reflectorArray(k)*changePhi; 
           end
       end
       
       
       function out = getReflectorArrays(obj)
          out = obj.reflectorArray; 
       end
       
       
   end
   
   %PRIVATES METHODS
   methods(Access = private)
       
     function ret = propagateSingleWave(obj,wave,initP,endP)
           ret = wave;
           increase = 1;
           if initP > endP
              increase = -1; 
           end
           for k = initP+1:increase:endP
               
              ret = ret*exp(obj.propagationModule.propagationArray(k).gamma);
          end
     end
       
   end
   
   
end