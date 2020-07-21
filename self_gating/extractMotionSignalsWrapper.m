function [Respiratory, Cardiac, respRange] = extractMotionSignalsWrapper(signal, opt)

method = opt.sigExt;

switch method
    
    case 'PCA'
        
        [Respiratory, Cardiac, respRange] = extractMotionSignals(signal, opt);
        
    case 'SSAFARI'
        
        [Respiratory, Cardiac] = extractMotionSignals_SSAFARI(signal, opt);
        respRange = [];
        
    case 'PT'
        
        [Respiratory, ~] = extractMotionSignals_PT(signal, opt);
        [~, Cardiac] = extractMotionSignals(signal, opt);
        respRange = [];
        
    case 'PTSola'
        [Respiratory, Cardiac] = extractMotionSignals_PTSola(signal, opt);
        respRange = [];
        
end

end

