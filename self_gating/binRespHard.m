function [rBins] = binRespHard(lengthDat, R, opt)

% Amplitude-based respiratory binning
% Range between 95th and 5th percentile to increase robustness to outliers
switch opt.sigExt
    case 'PCA'
%         dR = prctile(R,95) - prctile(R,5);
        dR = max(R) - min(R);
        
    case 'SSAFARI'
        dR = max(R) - min(R);
end

% Bin size
binWidth = dR / (opt.nRPhases); 

% Sinc interpolation of full data
rAmpDat = R;
% rAmpDat = (resampleSINC(rAmpDat',opt.sgI))';
rAmpDat = interp1(rAmpDat,1+(1/opt.sgI):(1/opt.sgI):length(rAmpDat)+1,'linear'); 
tmp = zeros(lengthDat,1);
tmp(1:opt.sgI-1) = ones(opt.sgI-1,1)*R(1); tmp(opt.sgI:end) = rAmpDat(1:end-opt.sgI+1);
rAmpDat = tmp;

% Respiratory bins
bin = zeros(1,opt.nRPhases);
for i = 1 : opt.nRPhases
%     bin(i) = prctile(R,95) - (i-1) * binWidth; 
    bin(i) = max(R) - (i-1)*binWidth;
end

% Start binning
for i = 1 : length(rAmpDat)
    for j = 1 : length(bin)
        if rAmpDat(i) <= bin(j)
            rBins(i) = j;
        elseif rAmpDat(i) > bin(1)
            rBins(i) = 1;         
        end
    end
end

% *************************************************************************
% plot bins
t = [opt.firstSample:1:opt.lastSample]*opt.TR;
figure;
for j = 1 : length(bin)
    
   plot([t(1),t(end)],[bin(j),bin(j)],'k'); hold on;
    
end

plot(t,R,'b'); hold off; %plot(t,R,'.k'); hold off;

figure;
plot(rAmpDat);
% *************************************************************************
