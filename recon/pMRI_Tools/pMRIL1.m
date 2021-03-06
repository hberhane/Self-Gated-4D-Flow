function [xhat,maps,GAMPout,param] = pMRIL1(kdata,samp,wvar,lambda,maps)

uniform_var = 0;
precision = 'single';
compute = 'mat';

% force the data to be downsampled
kdata = bsxfun(@times,kdata,permute(samp,[1,2,4,3]));

if nargin<5
    %% Estimate sensitivity maps
    weights_b = repmat(sum(samp,3),[1,1,size(kdata,3)]);
    weights_b(find(weights_b==0)) = Inf;

    % Normalize k_data by weights and take an inverse Fourier Transform
    time_av_b = ifft2_shift(sum(kdata,4)./weights_b);
    [x0, maps] = WalshCoilCombine(time_av_b,3);
    maps = repmat(maps,[1,1,1,size(samp,3)]);
    x0 = repmat(x0,[1,1,size(samp,3)]);
else
    %% Estimate sensitivity maps
    weights_b = repmat(sum(samp,3),[1,1,size(kdata,3)]);
    weights_b(find(weights_b==0)) = Inf;

    % Normalize k_data by weights and take an inverse Fourier Transform
    time_av_b = ifft2_shift(sum(kdata,4)./weights_b);
    [x0, ~] = WalshCoilCombine(time_av_b,3);
%     maps = repmat(maps,[1,1,1,size(samp,3)]);
    x0 = repmat(x0,[1,1,size(samp,3)]);
end

%% fftshift and downsample data
kdata = fftshift(fftshift(kdata,1),2);
samp = fftshift(fftshift(samp,1),2);
x0 = fftshift(fftshift(x0,1),2);
maps = fftshift(fftshift(maps,1),2);
kdata = downsample_data(kdata,samp);

% Normalize The columns of A to be unit norm
R = numel(samp)/length(find(samp ==1));
maps = maps*sqrt(R);
kdata = kdata*sqrt(R);
wvar = wvar*R;

use_gpu = 1;
if use_gpu
   kdata = gpuArray(single(kdata));
   samp = gpuArray(single(samp));
   x0 = gpuArray(single(x0));
   maps = gpuArray(single(maps));
   compute = 'gpu';
   precision = 'double';
end

%% Create Operators =======================================================
pMRI = pMRI_Op_2D_t(maps,samp,'uniform_var',uniform_var,'precision',precision,'compute',compute);

% Create nd-DWT Linear Transform
w_prop.wname = {'db1','db1','db1'};
w_prop.level = 1;
W = ndDWTLinTrans(w_prop,size(samp),'uniform_var', uniform_var,'compute',...
    compute,'precision',precision);

% Concatonate All The Linear Transform Operator Together
Op_b = LinTransConcat({pMRI;W},[1,1],precision,compute);

% Create Estimation classes
inputEst = NullEstimIn(0,1);
MeasEstimOut = CAwgnEstimOut(kdata,wvar,1);

%% Create a Vector of lambda values
% lambda = setLambda(size(samp),lambda);
k = 5;
lambda_band = zeros(8,1);
lambda_band(1) = 0.01*lambda;
lambda_band(2) = lambda;
lambda_band(3) = lambda;
lambda_band(4) = lambda;
lambda_band(5) = k*lambda;
lambda_band(6) = k*lambda;
lambda_band(7) = k*lambda;
lambda_band(8) = k*lambda;

lambda_vec = zeros([size(samp),8]);
for ind = 1:8
   lambda_vec(:,:,:,ind) = lambda_band(ind);
end

AnaEstimOut  = CplxLaplaceEstimOut(lambda_vec(:));
EstimOut = EstimOutConcat({MeasEstimOut;AnaEstimOut},[pMRI.M,W.M],precision,compute);

% Gamp Options
GAMPopt = GampOpt();
GAMPopt.nit = 100;
GAMPopt.adaptStep = true;
GAMPopt.stepIncr = 1.1;
GAMPopt.stepWindow = 0;
GAMPopt.step = 0.1;
GAMPopt.verbose = 1;
GAMPopt.legacyOut = 0;
GAMPopt.xhat0 = x0(:);


% Run GAMP
xhat = gampEst(inputEst,EstimOut,Op_b,GAMPopt);
    
GAMPout = xhat;
xhat = reshape(xhat.xhat,size(samp));
xhat = fftshift(fftshift(xhat,1),2);

xhat = gather(xhat);

end
