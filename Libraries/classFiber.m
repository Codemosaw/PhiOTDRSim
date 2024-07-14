classdef classFiber < handle
   properties
       n1               % Refraction index of the core
       n2               % Refraction index of the cladding
       coreRadio        % Radio of the core in meters
       fiberLength      % Length of the fiber in meters
       segmentLength    % Length of the segment in meters
       
       numberOfSegments % The total amount of segments in the fiber
   end
    
   methods
       
       function obj = classFiber(n1,n2,coreRadio,fiberLength,segmentLength)
       % Brief:
       %    Constructor for the fiber class, the class contain all the data
       %    regarding to the dimension and the means of the refraction
       %    index
       % Sintax:
       %    obj = classFiber(n1,n2,coreRadio,fiberLength,segmentLength)
       % inputs:
       %    * n1 - Refraction index of the core
       %    * n2 - Refraction index of the cladding 
       %    * coreRadio - Radio of the core in meters
       %    * fiberLength - Length of the fiber in meters
       %    * segmentLength - Length of the segment in meters
       % outpus:
       %    * obj - A classFiber instance
       % notes: 
       %    * fiberLength will be adjusted to be divisible for segmentLength
       %
       % methods:
       %    *   NONE
               if(n2>n1)
                  error("There can't be total internal reflection, n2 is greater than n1") 
               end
           
               obj.n1 = n1;
               obj.n2 = n2;
               obj.coreRadio = coreRadio;
               obj.segmentLength = segmentLength;
               
               obj.numberOfSegments = floor(fiberLength/segmentLength);
               
               obj.fiberLength = obj.numberOfSegments*segmentLength;
           end
   end
   
   %PRIVATES METHODS
   methods(Access = private)
       
     
       
   end
   
   
end