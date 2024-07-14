

fgrMng.newFigure();

for k = 1:length(range)
    
   PhaseK = unwrap(angle(saveMatrix(k,:)));
   
   plot(PhaseN) 
   hold on
   plot(PhaseK)
   xlim([420,920])
   grid on
   hold off
   
   pause(1);
   
end


