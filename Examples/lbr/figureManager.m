classdef figureManager < handle
    properties (Access = private)
       figureNumber 
    end
    
    methods 
        function obj = figureManager()
        % brief:
        %   This class is used only to manage figures in a isolated way
        % notes:
        %   This function has nothing to do whit the PhiOTDR simulator
            obj.figureNumber = 1;
        end
        
        function fgr = newFigure(obj)
           fgr = figure(obj.figureNumber);
           clf(fgr)
           obj.addCounter();
        end
        
        function obj = addCounter(obj)
            obj.figureNumber = obj.figureNumber + 1; 
        end
        
        function obj = resetCounter(obj)
           obj.figureNumber = 1; 
        end
    end
    
end