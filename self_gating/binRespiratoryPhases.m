function [rBins, weights] = binRespiratoryPhases(lengthDat, R, opt)
% =========================================================================
% This is a test function which creates 1 respiratory bin surrounding
% end-expiration and assigns Gaussian weights to each point within the bin
% =========================================================================

% Find the center of the appropriate bin (end-expiration)
tmp = R;
counter = 0;
while(1)
    counter = counter + 1;
    nL = length(find(tmp < mean(tmp)));
    nG = length(tmp) - nL;
    
    if counter == 1
        if nL >= nG
            direction = 'ascend';
            tmp(tmp > mean(tmp)) = [];
        elseif nL < nG
            direction = 'descend';
            tmp(tmp < mean(tmp)) = [];
        end
        
    elseif counter > 1
        if nL > nG
            tmp(tmp > mean(tmp)) = [];
        elseif nL < nG
            tmp(tmp < mean(tmp)) = [];
        else
            % Center of desired (most populated) respiratory phase
            cPh = mean(tmp);
            break;
        end
    end
end

% Find upper boundary of the bin
[sortR, indSort] = sort(R, direction);
cdf = zeros(length(R),1);
for i = 1 : length(R)
    cdf(i) = mean(sortR(1:i) - cPh);
end

% Bin edges
coInd = find(abs(cdf) == min(abs(cdf)));
% binLo = min(sortR(1), sortR(coInd));
% binHi = max(sortR(1), sortR(coInd));

% Resample to match data
rAmpDat = R;
% rAmpDat = resample(rAmpDat,N,D);
rAmpDat = (resampleSINC(rAmpDat',opt.sgI))';
tmp = zeros(lengthDat,1);
tmp(1:opt.sgI-1) = ones(opt.sgI-1,1)*R(1); tmp(opt.sgI:end) = rAmpDat(1:end-opt.sgI+1);
rAmpDat = tmp;

% Trim self gating signal out of resampled amplitude data
% rAmpDat(opt.sgI:opt.sgI:end) = [];
% mu = cPh;
% sigma = std(R);
% Z = (1/(sigma*sqrt(2*pi)))*exp(-((R-mu).^4)/(2*sigma^4));
% FWHM = 2*nthroot(-2*(sigma^4)*log(0.5),4);

FWHM = std(R);
sigma = (FWHM/2)*nthroot((-1)/(2*log(0.5)),4);
mu = cPh;
Z = (1/(sigma*sqrt(2*pi)))*exp(-((R-mu).^4)/(2*sigma^4));
Z_sort = (1/(sigma*sqrt(2*pi)))*exp(-((sortR-mu).^4)/(2*sigma^4));

binLo = cPh - (sigma/2);
binHi = cPh + (sigma/2);

% =========================================================================
% Plot bins
figure;
t = [1:1:length(R)] * opt.TR;
plot([1,t(end)],[binLo,binLo],'--k'); hold on;
plot([1,t(end)],[binHi,binHi],'--k'); hold on;
plot(t,R,'b'); hold on; %plot(t,R,'.k'); hold off;
% plot([1,t(end)],[cPh,cPh],'k');
% plot([1,t(end)],[cPh+(FWHM/2),cPh+(FWHM/2)],'--r'); hold on;
% plot([1,t(end)],[cPh-(FWHM/2),cPh-(FWHM/2)],'--r'); hold on;
% plot([1,t(end)],[cPh+(sigma/2),cPh+(sigma/2)],'-r'); hold on;
% plot([1,t(end)],[cPh-(sigma/2),cPh-(sigma/2)],'-r'); hold on;

% =========================================================================

% Bin the k-space lines
rBins = zeros(lengthDat, 1);
for i = 1 : length(rAmpDat)
    if rAmpDat(i) >= binLo && rAmpDat(i) <= binHi
        rBins(i) = 1;
    end
end

% Assign Guassian weights within bin
tmp = rAmpDat(logical(rBins));
pd = fitdist(tmp,'Normal');
weights = pdf(pd, rAmpDat);
weights = weights / sum(weights);

rBins(rBins==0) = 2;


