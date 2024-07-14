
classdef classPropagationSegment < handle
   properties
       attenuation;            
       phi;             %phase change in this segment  
       gamma;
       n;
       beta;
       
   end
    
   methods
       % brief:
       %    Constructor for the propagation segment class
       %
       % sintax:
       %     obj = classPropagationSegment(attenuation, phi,n,beta)
       % Description:
       %    Create a segment of propagation, used by the classPropagation class
       %
       % inputs:
       %    attenuation - integration of the attenuation between the whole segment
       %    phi - integration of the phase change between the whole segment
       %    n - relexion index in the whole segment (the average one)
       %    beta - value of beta in the whole segment (the average one)
       %
       % note:
       %    * this class is intended to be used by the classPropagator, the user
       %      SHOULD NOT use it directly

       function obj = classPropagationSegment(attenuation, phi,n,beta)
               obj.attenuation = attenuation;
               obj.phi = phi;
               obj.gamma = attenuation + 1i*phi;
               
               obj.n = n;
               obj.beta = beta;
        end
   end
   
   %PRIVATES METHODS
   methods(Access = private)
       
     
       
   end
   
   
end