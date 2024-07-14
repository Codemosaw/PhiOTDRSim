classdef classPerturbator < handle
   properties(Access = private)
    gamma_T % The thermo-Optic Coefficient 
    eta_T % The thermal-Expansion Coefficient
    
    u_s % The Poisson Ratio
    p12 % The tenso-optic Coefficient 1
    p11 % The tenso-optic Coefficient 2
    
    epsilon_e % The tension coefficient used in the calculations of phase variations
    
    propagador % The propagation instance used

   end
    
   methods

       function obj = classPerturbator(gamma_T,eta_T,u_s,p12,p11,propagador)
              % brief: 
              %     Create a instance of the perturbator class to set the perturbations
              % sintax:
              %     obj = classPerturbator(gamma_T,eta_T,u_s,p12,p11,propagador)
              % Description:
              %     Create the perturbation class
              % inputs:
              %     * gamma_T - thermo-Optic Coefficient
              %     * eta_T - thermal-Expansion Coefficient
              %     * u_s - Poisson Ratio
              %     * p12 - tenso-optic Coefficient 1
              %     * p11 - tenso-optic Coefficient 2
              %     * propagador - A classPropagation class instance that 
              %                    will be linked to this classPerturbator
              %                    instance
              % outputs:
              %     * obj - A classPerturbator instance
              %
              %
              % Methods (Use command help whit each method to learn more) :
              %     * temperatureChange(obj,deltaTemperature, initialPoint, finalPoint)
              %     * strainChange(obj,strain, initialPoint, finalPoint)
              
              obj.eta_T = eta_T;
              obj.gamma_T = gamma_T;
              obj.u_s =u_s;
              obj.p12 = p12;
              obj.p11 = p11;
              obj.propagador = propagador;
              obj.epsilon_e = -(1/2)*(propagador.n^2)*( (1-obj.u_s)*obj.p12-obj.u_s*obj.p11);
       end
       
       
       
       function nT = temperatureChange(obj,deltaTemperature, initialPoint, finalPoint)
       % brief:
       %    Add Temperature perturbation to the linked classPropgation
       %    class instance
       % sintax:
       %    nT = temperatureChange(deltaTemperature, initialPoint, finalPoint)
       % inputs:
       %    * deltaTemperature - the variation of temperature, in kelvin
       %                         degress
       %    * initialPoint - initial distance point in which the variation
       %                     will be added, in meters
       %    * finalPoint - final distance point in which the variation will
       %                   be added, in meters
       % returns:
       %    * nT - The index refrection variations equivalent to the
       %           temperature variation
       % notes:
       %    * Thermo-optic coefficient referencial value 5*10^-6
       %    * Thermal expansion coefficient referencial value 5*10^-7
       
       iPoint = initialPoint/obj.propagador.fiber.segmentLength;
       fPoint = finalPoint/obj.propagador.fiber.segmentLength;
           nT = deltaTemperature*(obj.gamma_T + obj.propagador.n*obj.eta_T);
           
           obj.propagador.changeRefractiveIndex(nT,iPoint,fPoint);
           
       end
       
       

       function ne = strainChange(obj,strain, initialPoint, finalPoint)
       % brief:
       %    Add strain perturbation to the linked classPropgation
       %    class instance
       % sintax:
       %    ne = strainChange(strain, initialPoint, finalPoint)
       % inputs:
       %    * strain - the variation of strain
       %    * initialPoint - initial distance point in which the variation
       %                     will be added, in meters
       %    * finalPoint - final distance point in which the variation will
       %                   be added, in meters
       % returns:
       %    * ne - The index refrection variations equivalent to the
       %           strain variation
           ne = strain*(obj.epsilon_e+1)*obj.propagador.n;
           
           iPoint = initialPoint/obj.propagador.fiber.segmentLength;
           fPoint = finalPoint/obj.propagador.fiber.segmentLength;

           obj.propagador.changeRefractiveIndex(ne,iPoint,fPoint);
       end

   end
   
   %PRIVATES METHODS
   methods(Access = private)
       
     
       
   end
   
   
end