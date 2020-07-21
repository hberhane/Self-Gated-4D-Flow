function [Respiratory, Cardiac, respRange] = extractMotionSignals(SG, opt)
%% OPTIONS
if ~isfield(opt,'nCoils');  opt.nCoils  = 12;           end     % Coil elements used for PCA (not used right now!)
if ~isfield(opt,'rLPF');    opt.rLPF    = [0.5];    end     % FIR low pass filter to select respiratory signal (Hz)
if ~isfield(opt,'cBPF');    opt.cBPF    = [0.5 3];      end     % FIR band pass filter to select cardiac signal (Hz)
if ~isfield(opt,'rcLPF');   opt.rcLPF   = [3];     end     % FIR low/band pass filter to select both respiratory and cardiac signal (Hz)
if ~isfield(opt,'nCHarm');  opt.nCHarm  = 2;            end     % # of cardiac harmonics to use (not implemented at the moment!)
if ~isfield(opt,'fOrd');    opt.fOrd    = 1001;         end     % FIR filter order (length+1), introudces transients so lower=better assuming suffiecient cuttoffs
if ~isfield(opt,'fs');      opt.fs      = 20;           end     % Self-gating sampling frequency in Hz

%% PREPROCESSING
for j = 1:size(SG,2) 
    for k = 1 : size(SG,3)        
        tmp = squeeze(SG(:,j,k))';
        % Reconstruct 1D projections of magnitude 
        A(:,j,k) = double(abs(ifftshift(ifft(fftshift(tmp),[]))));
%         A(:,j,k) = double((ifftshift(ifft(fftshift(tmp),[]))));
    end
end

% A = cat(2, real(A), imag(A));

% % Little piece code for pilot tone
% % ********************************
% ind = round(size(A,1)/4);
% A_pt(1:ind,:,:) = A(1:ind,:,:);
% A_pt(ind+1:2*ind,:,:) = A(3*ind+1:size(A,1),:,:);
% A_pt_sum = sum(abs(A_pt),3);
% A_pt_SOS = sum(A_pt_sum.^2,2);
% 
% pt_ind = find(A_pt_SOS == max(A_pt_SOS));
% 
% pt_sig = A_pt(pt_ind,:,:);
% 
% A = pt_sig;
% *********************************
% for j = 1 : size(SG, 2)
%     for k = 1 : size(SG, 3)
%     
%         A(1,j,k) = SG(floor((size(SG,1)/2))+1, j, k);
%         
%     end
% end
% A = abs(A);

A_orig = A;
% Concatenate coil elements to readout dimension
A = reshape(A,size(A,1)*size(A,2),size(A,3));

% Tranpose A for PCA -> [Time]x[Readout]
A = A'; 

B = A;
% **********
% blah = A(:, 1800:1950)';
% blah = blah(55:80,:);

% Spectral analysis
% *****************
% n = size(A,1);
% df = 1/(n*opt.TR);
% f = [-n/2:1:n/2-1]*df;
% for i = 1 : 1000
%     y(:,i) = abs(fftshift(fft(A(:,i))));
% end
% ymean = mean(y, 2);
% plot(f,y);
% ********************

% Display matrix A
figure; imagesc(A);

% Subtract the mean from each column of A
A = A - repmat(mean(A,1),[size(A,1),1]);

% FIR filtering over temporal dimension to select respiratory/cardiac signals
rFilt = fir1(opt.fOrd-1,opt.rLPF/(opt.fs/2),hanning(opt.fOrd));     % Respiratory
cFilt = fir1(opt.fOrd-1,opt.cBPF/(opt.fs/2),hanning(opt.fOrd));     % Cardiac
% rcFilt = fir1(fOrd-1,rcLPF/(fs/2),hanning(fOrd));   % Both
for j =  1 : size(A,2)
    % Select respiratory signal
    rA(:,j) = conv(A(:,j),rFilt,'same');
    % Select cardiac signal
    cA(:,j) = conv(A(:,j),cFilt,'same');
%     % Select both  
%     rcA(:,j) = conv(A(:,j),rcFilt,'same');
    % 'same' flag corrects for group delay caused by FIR filtering
end

%% AUTOMATIC COIL SELECT (NOT IMPLEMENTED AT THE MOMENT!)
% 
% A_orig2 = A'; A_orig2 = reshape(A_orig2,[size(A_orig,1),size(A_orig,2),size(A_orig,3)]);
% channelInd = autoSelectChannel(A_orig, opt); 
% 
% % auto channel select for best respiratory signal
% % colNorm = sqrt(sum(rA.^2,1));
% % cInd = find(colNorm == max(colNorm));
% % 
% % low = cInd - 64;
% % hi = cInd + 64;

%% EXTRACT RESPIRATORY AND CARDIAC MOTION USING PCA
% Principle component analysis 
% Respiratory signal
Respiratory = pcaSG(rA,1);

% Cardiac signal
Cardiac = pcaSG(cA,1);

%% ESTIMATE RESPIRATORY RANGE FOR WEIGHTING

% Find coil with highest correlation with respiratory signal
for i = 1 : size(A, 2)
    tmp = corrcoef(Respiratory', rA(:,i));
    Ccoef(i) = tmp(2,1);
end

% tmp = Ccoef;
% tmp(tmp<=0.0) = 0;
% tmp = logical(tmp);
% start1 = strfind([0,tmp==1],[0 1]);
% end1 = strfind([tmp==1,0],[1 0]);
% width1 = end1 - start1;
% width1(width1==0) = [];
% width1(width1==1) = [];
% width1(width1==2) = [];
% width1(width1==3) = [];
% width1(width1==4) = [];
% 
% tmp = -Ccoef;
% tmp(tmp<=0.0) = 0;
% tmp = logical(tmp);
% start1 = strfind([0,tmp==1],[0 1]);
% end1 = strfind([tmp==1,0],[1 0]);
% width2 = end1 - start1;
% width2(width2==0) = [];
% width2(width2==1) = [];
% width2(width2==2) = [];
% width2(width2==3) = [];
% width2(width2==4) = [];
% 
% width = cat(2,width1,width2);
% respRange = mean(width);


% tmp = Ccoef;
% for i = 1 : 12
%     pkInd(i) = find(abs(tmp)==max(abs(tmp)));
%     coilNum(i) = ceil(pkInd(i) / size(SG,1));
%     tmp(size(SG,1)*(coilNum(i)-1)+1:size(SG,1)*coilNum(i)) = 0;
%     coil(:,:,i) = squeeze(A_orig(size(A_orig,1)/4+1:end-size(A_orig,1)/4,coilNum,:));
% end

% for i = 1 : length(pkInd)
%     
%     if Ccoef(pkInd(i)) >= 0
%         tmp = Ccoef;
%     else
%         tmp = -Ccoef;
%     end
%     
%     tmp(tmp <= 0.5) = 0;
%     Ccoef2 = tmp;
%     
%     cInd = pkInd(i);
%     lowInd = cInd;
%     while 1
%         if Ccoef2(lowInd) <=0 || lowInd == size(SG,1)*(coilNum(i)-1)+1
%             break
%         end
%         lowInd = lowInd - 1;
%     end
%     hiInd = cInd;
%     while 1
%         if Ccoef2(hiInd) <= 0 || hiInd == size(SG,1)*coilNum(i)
%             break
%         end
%         hiInd = hiInd + 1;
%     end
%     
%     respRange(i) = (hiInd-1) - (lowInd+1) + 1;
%     
% end
% 
% zscore = (respRange-mean(respRange))/std(respRange);
% respRange(abs(zscore)>2) = [];
% 
% respRange_mean = mean(respRange);


    
pkInd = find(abs(Ccoef)==max(abs(Ccoef)));
coilNum = ceil(pkInd / size(SG,1));
coil = squeeze(A_orig(size(A_orig,1)/4+1:end-size(A_orig,1)/4,coilNum,:));
% apply respiratory filter
for i = 1 : size(coil, 1)
    fCoil(i,:) = conv(coil(i,:),rFilt,'same');
end
coil = fCoil;

% determine which coef. sign dominates and threshold the rest
for i = 1 : size(coil, 1)
    tmp = corrcoef(Respiratory', coil(i,:));
    Ccoef2(i) = tmp(2,1);
end
% pkInd = find(abs(Ccoef2)==max(abs(Ccoef2)));
% if Ccoef2(pkInd) > 0
%     coil(Ccoef2<0,:) = 0;
% else
%     coil = -coil; Ccoef2 = -Ccoef2;
%     coil(Ccoef2<0,:) = 0;
% end
if sum(Ccoef2(Ccoef2>0)) >= sum(Ccoef2(Ccoef2<0))
    coil(Ccoef2<0,:) = 0;
else
    coil(Ccoef2>0,:) = 0;
    coil = -coil; Ccoef2 = -Ccoef2;
end

%extract most pertinent band
% width of band is respiratory range in pixels
cInd = find(Ccoef2 == max(Ccoef2));
lowInd = cInd;
while 1
    if Ccoef2(lowInd) <=0 || lowInd == 1
        break
    end
    lowInd = lowInd - 1;
end
hiInd = cInd;
while 1
    if Ccoef2(hiInd) <= 0 || hiInd == length(Ccoef2)
        break
    end
    hiInd = hiInd + 1;
end

coil([1:lowInd,hiInd:end],:) = 0;
respRange = (hiInd-1) - (lowInd+1);

figure;
subplot(221);   imagesc(B); hold on
                plot([size(SG,1)*(coilNum-1)+1,size(SG,1)*(coilNum-1)+1],[1,size(A,1)],'-k','LineWidth',2);
                plot([size(SG,1)*coilNum,size(SG,1)*coilNum],[1,size(A,1)],'-k','LineWidth',2);
                title(['Casorati']);
subplot(222);   imagesc(fCoil);
                title(['Selected Coil: ', num2str(coilNum)]);
subplot(223);   plot(Ccoef2);
                title('Correlation coefficients');
subplot(224);   imagesc(coil);
                title(['Range = ', num2str(respRange), ' pixels']);


% for i = 1 : size(coil, 2)
%     cSum = cumsum(coil(:,i));
%     cSum = cSum - min(cSum); cSum = cSum / max(cSum);
%     x1 = find(cSum == max(cSum(cSum < 0.5)));
%     y1 = cSum(x1); y2 = cSum(x1+1);
%     
%     centerOfMass(:,i) = x1 + ((0.5 - y1) / (y2 - y1));
% end
% respRange = max(centerOfMass) - min(centerOfMass);

end

