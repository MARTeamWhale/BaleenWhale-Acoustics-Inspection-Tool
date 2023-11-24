% function for saving figure based on pixel dimensions. Because MATLAB
% unecessarily complicates things.
% Currently only supports PNG format.
% Sep 2016

function saveFigInPixels(hfig,figPath,pixelDim)
    
    % 1) get relevant figure properties
    oldPU = hfig.PaperUnits;
    oldPP = hfig.PaperPosition;
    oldPPM = hfig.PaperPositionMode;
    
    % 2) change figure properties
    hfig.PaperUnits = 'inches';
    hfig.PaperPosition = [0,0,pixelDim/150];
    hfig.PaperPositionMode = 'manual';
    
    % 3) print figure
    print(hfig,figPath,'-dpng','-r150')
    
    % 4) return properties back to old
    hfig.PaperUnits = oldPU;
    hfig.PaperPosition = oldPP;
    hfig.PaperPositionMode = oldPPM;
end