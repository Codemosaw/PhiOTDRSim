% author: Marcos/Sasmosaw Rojas Mardones
classdef classTransmitter < handle
   properties
       lambda               % wavelength
       field                % incident length, its a complex number
       widthSeconds               % width of the pulse in seconds

   end
    
   methods
       function obj = classTransmitter(lambda, field, widthS)
       % brief:
       %    Constructor for the tranmitter class
       % sintax:
       %    obj = classTransmitter(lambda, field, widthS)
       % inputs:
       %	* lambda - Wavelength in the void of the input, in m (units
       %               are adjusted inside)
       %    * field - complex number representing the phase and amplitud
       %              of the electric field
       %	* widthS - Width of the pulse in seconds
       % notes:
       %    * The widthS parameter will be aproximatted for simplicity,
       %      precission will be lost
       %
       % methods:
       %    * NONE
               obj.lambda = lambda;
               obj.field = field;
               obj.widthSeconds = widthS;
           end
   end
   
   %PRIVATES METHODS
   methods(Access = private)
       
     
       
   end
   
   
end