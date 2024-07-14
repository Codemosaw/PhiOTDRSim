
classdef classPropagation < handle
   properties
       alfa
       beta
       b %normalized beta
       n %effective refrective index
       V % normalized frecuency
       k0 %wave number
       Vp %propagation
       
       betaFactor %to adjust
       alfaFactor %to adjust
       
       gamma
       
       fiber
       
       emitter 
       widthSegments  %width of the pulse in number of segments
       
       
       propagationArray
   end
    
   methods
       
       function obj = classPropagation(fiber, emitter, alfaFactor, betaFactor)
       % brief:
       %    Constructor for the fiber class
       % sintax:
       %    obj = classPropagation(fiber, emitter, alfaFactor, betaFactor)
       % inputs:
       %    * fiber - Must be a classFiber instance
       %    * emitter - Must be a classTransmitter instance
       %    * alfaFactor - Its a simple gain for the attenuation, used to
       %                 adjust the losses of the fiber
       %    * betaFactor - Its a simple gain for the beta parameter,
       %    normally it should be settd as zero
       % outpus:
       %    * obj - A classPropagation instance
       % notes:
       %    * A ideal fiber can be simulated whit alfaFactor equal to 0
       %
       % methods:
       %    * propagateWave(initP,endP)
       %    * propagateWaveSaveArray(initP,endP)
       %    * changeRefractiveIndex(deltaN,initialSegment,finalSegment)
       %    * resetPropagationArrays()
       
            obj.fiber = fiber;
            obj.emitter = emitter;
            
            %Extrayendo datos fundamentales de la fibra
            obj.k0 = 2*pi/emitter.lambda;
            obj.V = obj.k0*obj.fiber.coreRadio*sqrt(fiber.n1^2-fiber.n2^2);
            
            if obj.V > 2.405 
                warning("Out of monomode domain,aproximatoin may not be accuerated")
            end
            
            %aproximación basica de b en función de V
            V = obj.V;
            obj.b = -0.0574*V.^4 + 0.2544*V.^3 - 0.2309*V.^2 + 0.1001*V - 0.0036;
            b = obj.b;
            obj.n = fiber.n2 - b*(fiber.n2 - fiber.n1);
            
            obj.beta = betaFactor*obj.k0*sqrt( b*(obj.fiber.n1^2 - obj.fiber.n2^2) + obj.fiber.n2^2 );
            obj.betaFactor = betaFactor;
            
            obj.alfa = alfaFactor*(1.7/(4.343*1000))*(0.85/(emitter.lambda/10^(-6)))^4;
            obj.alfaFactor  = alfaFactor;
            
            obj.gamma = obj.alfa + 1i*obj.beta;
            
            obj.Vp = 3*10^8/obj.n;
            
            %array for each segment
            obj.propagationArray = classPropagationSegment.empty(0,fiber.numberOfSegments);
            
            segmentLength = obj.fiber.segmentLength;
            for n=1:obj.fiber.numberOfSegments
                obj.propagationArray(n) = classPropagationSegment(-obj.alfa*segmentLength,-obj.beta*segmentLength,obj.n,obj.beta);
            end
            
            
            %width of the pulse in segments
            segments = obj.emitter.widthSeconds*obj.Vp;
            widthSegments = floor(segments/obj.fiber.segmentLength);
            
            if widthSegments < 2
               widthSegments = 2;
               warning("pulse is to short, it will be set to 2, results may not be accuearate")
            end
            
            obj.widthSegments = widthSegments;
            
       end
       
       
       
       
       function ret = propagateWave(obj,initP,endP)
       % brief:
       %    This functions propagate a wave from a given initial point,
       %    to a final point
       % sintax:
       %    ret = propagateWave(initP,endP)
       % inputs:
       %    * initP - Initial point (segment) in which the wave will be
       %              propagated
       %    * endP - Final point (segment) in which the wave will be
       %             propagated
       % notes:
       %    * All the segments points are measured from the border of the
       %      segment, for example a fiber whit 3 segment will have points from
       %      zero to the number of segments of the fiber
           wave = obj.emitter.field;
           ret = wave;
           increase = 1;
           if initP > endP
              increase = -1; 
           end
           for k = initP+1:increase:endP+1
              ret = ret*exp(obj.propagationArray(k).gamma);
          end
       end
       
       
       
       
       function [d,t,E] = propagateWaveSaveArray(obj,initP,endP)
       % brief:
       %    This functions propagate a wave from a given initial point,
       %    to a final point
       % sintax:
       %    [d,t,E] = propagateWaveSaveArray(initP,endP)
       % inputs:
       %    *initP - Initial point (segment) in which the wave will be
       %             propagated
       %	* endP - Final point (segment) in which the wave will be
       %             propagated
       % outpus:
       %    * d Distance vector
       %	* t Time vector
       %	* E Complex field vector
           wave = obj.emitter.field;
           increase = 1;
           if initP > endP
              increase = -1; 
           end
           E = zeros(1,length(initP:increase:endP));
           
           E(1) = wave;
           for k = 1:length(E)-1
              E(k+1) = E(k)*exp(obj.propagationArray(k).gamma);
           end
           
           d = (initP:increase:endP)*obj.fiber.segmentLength;
           t = d/obj.Vp;
           
       end
       
       
       
       function obj = changeRefractiveIndex(obj,deltaN,initialSegment,finalSegment)
       % brief:
       %    change the effecitve index por the fiber
       % sintax:
       %    changeRefractiveIndex(deltaN,initialSegment,finalSegment)
       % inputs:
       %	* deltaN - Change in the refrective index
       %	* initalSegment - Is the segment number in which the
       %                      perturbation begins
       %	* finalSegment - Is the segment number in which the
       %                     perturbation ends
       % outpus: NONE
          for k = initialSegment:finalSegment
            
            deltaB = deltaN*obj.k0;
            deltaPhi = deltaB*obj.fiber.segmentLength;
            
            obj.propagationArray(k).n = obj.propagationArray(k).n + deltaN;
            obj.propagationArray(k).gamma = obj.propagationArray(k).gamma - 1i*deltaPhi;
            obj.propagationArray(k).phi = obj.propagationArray(k).phi - deltaPhi;
            obj.propagationArray(k).beta = obj.propagationArray(k).beta + deltaB;
            
          end
          
       end
       
       
       function obj = resetPropagationArrays(obj)
       % brief:
       %    Reset all the propagatio array used to restar the propagation
       %    instance to its initial state
       % sintax:
       %    resetPropagationArrays()
       % inputs: NONE
       % outpus: NONE
            
            obj.propagationArray = classPropagationSegment.empty(0,obj.fiber.numberOfSegments);
            
            segmentLength = obj.fiber.segmentLength;
            for k=1:obj.fiber.numberOfSegments
                obj.propagationArray(k) = classPropagationSegment(-obj.alfa*segmentLength,-obj.beta*segmentLength,obj.n,obj.beta);
            end
            
            
       end
       
   end
   
   %PRIVATES METHODS
   methods(Access = private)
       
     
       
   end
   
   
end