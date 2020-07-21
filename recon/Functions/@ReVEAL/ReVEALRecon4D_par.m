function ReVEALRecon4D_par(obj)

% Allow parallel computing on 2 gpu devices simultaneously

%% fftshifts
applyFftshift(obj)

% obj.options.ReVEALOpts.gamma = fftshift(fftshift(obj.options.ReVEALOpts.gamma,1),2);
% obj.options.ReVEALOpts.gamma = obj.options.ReVEALOpts.gamma(:);

%% Normalize The columns of A to be unit norm
if ~isempty(obj.data.weightsB)
    R = numel(obj.data.sampB)/sum(obj.data.weightsB(:).^2);
else  
    R = numel(obj.data.sampB)/length(find(obj.data.sampB ==1));
end
obj.data.sensMaps = obj.data.sensMaps*sqrt(R);
obj.data.Yb = obj.data.Yb*sqrt(R);
obj.data.Yx = obj.data.Yx*sqrt(R);
obj.data.Yy = obj.data.Yy*sqrt(R);
obj.data.Yz = obj.data.Yz*sqrt(R);
obj.options.ReVEALOpts.wvar = obj.options.ReVEALOpts.wvar*R;
obj.options.ReVEALOpts.sigmaSq = obj.options.ReVEALOpts.sigmaSq*R;
fprintf(sprintf('R = %s\n',num2str(R)))

% Check if multiple sensitivity maps are used from espirit
if size(obj.data.sensMaps, 5) == 2
    % Repeat data along time dimension
    % Data
    obj.data.Yb = repmat(obj.data.Yb,[1,1,1,1,2]);
    obj.data.Yx = repmat(obj.data.Yx,[1,1,1,1,2]);
    obj.data.Yy = repmat(obj.data.Yy,[1,1,1,1,2]);
    obj.data.Yz = repmat(obj.data.Yz,[1,1,1,1,2]);
    % Initialization
    obj.data.x0 = repmat(obj.data.x0,[1,1,1,2]);
    % Samples
    obj.data.sampB = repmat(obj.data.sampB,[1,1,1,2]);
    obj.data.sampX = repmat(obj.data.sampX,[1,1,1,2]);
    obj.data.sampY = repmat(obj.data.sampY,[1,1,1,2]);
    obj.data.sampZ = repmat(obj.data.sampZ,[1,1,1,2]);
    % Weights
    if ~isempty(obj.data.weightsB)
    obj.data.weightsB = repmat(obj.data.weightsB,[1,1,1,2]);
    obj.data.weightsX = repmat(obj.data.weightsX,[1,1,1,2]);
    obj.data.weightsY = repmat(obj.data.weightsY,[1,1,1,2]);
    obj.data.weightsZ = repmat(obj.data.weightsZ,[1,1,1,2]);
    end
end

% ****************************
% Weight the data
if ~isempty(obj.data.weightsB)
   obj.data.Yb = bsxfun(@times,permute(obj.data.weightsB,[1,2,3,5,4]) ,obj.data.Yb);
   obj.data.Yx = bsxfun(@times,permute(obj.data.weightsX,[1,2,3,5,4]) ,obj.data.Yx);
   obj.data.Yy = bsxfun(@times,permute(obj.data.weightsY,[1,2,3,5,4]) ,obj.data.Yy);
   obj.data.Yz = bsxfun(@times,permute(obj.data.weightsZ,[1,2,3,5,4]) ,obj.data.Yz);   
end
% ***************************




% Set initializations
obj.options.GAMPOptB.xhat0 = obj.data.x0(:);
obj.options.GAMPOptX.xhat0 = obj.data.x0(:);
obj.options.GAMPOptY.xhat0 = obj.data.x0(:);
obj.options.GAMPOptZ.xhat0 = obj.data.x0(:);

% Downsample Data
dataB = downsample_data(obj.data.Yb,obj.data.sampB);
dataX = downsample_data(obj.data.Yx,obj.data.sampX);
dataY = downsample_data(obj.data.Yy,obj.data.sampY);
dataZ = downsample_data(obj.data.Yz,obj.data.sampZ);


% Create Input Estimation Class
pMRI_b = pMRI_Op_3D_t(obj.data.sensMaps,obj.data.sampB,'uniform_var',...
         obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute...
         ,'precision',obj.options.ReVEALOpts.precision,'weights',obj.data.weightsB);
pMRI_x = pMRI_Op_3D_t(bsxfun(@times,obj.data.sensMaps,obj.data.maxwellCorrX),obj.data.sampX,'uniform_var',...
         obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute...
         ,'precision',obj.options.ReVEALOpts.precision,'weights',obj.data.weightsX);
pMRI_y = pMRI_Op_3D_t(bsxfun(@times,obj.data.sensMaps,obj.data.maxwellCorrY),obj.data.sampY,'uniform_var',...
         obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute...
         ,'precision',obj.options.ReVEALOpts.precision,'weights',obj.data.weightsY);
pMRI_z = pMRI_Op_3D_t(bsxfun(@times,obj.data.sensMaps,obj.data.maxwellCorrZ),obj.data.sampZ,'uniform_var',...
         obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute...
         ,'precision',obj.options.ReVEALOpts.precision,'weights',obj.data.weightsZ);

% Create nd-DWT Linear Transform
W = ndDWTLinTrans4D(obj.options.SparseTrans,size(obj.data.sampB),'uniform_var',...
    obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute,...
    'precision',obj.options.ReVEALOpts.precision);

% Concatonate All The Linear Transform Operator Together
Op_b = LinTransConcat({pMRI_b;W},[1,1],obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
Op_x = LinTransConcat({pMRI_x;W},[1,1],obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
Op_y = LinTransConcat({pMRI_y;W},[1,1],obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
Op_z = LinTransConcat({pMRI_z;W},[1,1],obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);

% Create Output Estimation classes
MeasEstimOut_b = CAwgnEstimOut(dataB,obj.options.ReVEALOpts.wvar,obj.options.ReVEALOpts.MAP);
MeasEstimOut_x = CAwgnEstimOut(dataX,obj.options.ReVEALOpts.wvar,obj.options.ReVEALOpts.MAP);
MeasEstimOut_y = CAwgnEstimOut(dataY,obj.options.ReVEALOpts.wvar,obj.options.ReVEALOpts.MAP);
MeasEstimOut_z = CAwgnEstimOut(dataZ,obj.options.ReVEALOpts.wvar,obj.options.ReVEALOpts.MAP);

% Set the Regularization per band
lambda = setLambda(size(obj.data.sampB),obj.options.ReVEALOpts.lambda0);
AnaEstimOut1 = CplxLaplaceEstimOut(lambda);

% Create Output Estimation class
EstimOut_b = EstimOutConcat({MeasEstimOut_b;AnaEstimOut1},[pMRI_b.M,W.M],...
             obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
EstimOut_x = EstimOutConcat({MeasEstimOut_x;AnaEstimOut1},[pMRI_x.M,W.M],...
             obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
EstimOut_y = EstimOutConcat({MeasEstimOut_y;AnaEstimOut1},[pMRI_y.M,W.M],...
             obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
EstimOut_z = EstimOutConcat({MeasEstimOut_z;AnaEstimOut1},[pMRI_z.M,W.M],...
             obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);

% Create Input Estimation class
EstimIn = NullEstimIn(0,1);

% Clear some memory up
clear meanB meanX meanY meanZ varB varX varY varZ pvarb pvarx pvary pvarz...
      MeasEstimOut_b MeasEstimOut_x MeasEstimOut_y MeasEstimOut_z...
      AnaEstimOut1 pMRI_b pMRI_x pMRI_y pMRI_z W dataB dataX dataY dataZ ...
      lambda;

%% =========================================================================
% Initial Reconstruction 
tic


% ********* TEST OF PARALLEL COMPUTING ******************

S{1}.EstimIn = EstimIn;
S{1}.EstimOut = EstimOut_b;
S{1}.Op = Op_b;
S{1}.opt = obj.options.GAMPOptB;
S{1}.compute = 'cpu';

S{2}.EstimIn = EstimIn;
S{2}.EstimOut = EstimOut_x;
S{2}.Op = Op_x;
S{2}.opt = obj.options.GAMPOptX;
S{2}.compute = 'cpu';

pool = parpool(2);


% EstimI{1} = EstimIn;
% EstimI{2} = EstimIn;
% 
% EstimOut{1} = EstimOut_b;
% EstimOut{2} = EstimOut_x;
% 
% Op{1} = Op_b;
% Op{2} = Op_x;
% 
% opt{1} = obj.options.GAMPOptB;
% opt{2} = obj.options.GAMPOptX;

spmd
    if labindex == 1
        gpuDevice(1);
        % Reconstruct background image
        fprintf('\npB...\n')
        Out = reconGAMP(S{1});
    elseif labindex == 2
        gpuDevice(2);
        % Reconstruct x encoded image
        fprintf('\npX...\n')
        Out = reconGAMP(S{2});
    end
end


% Reconstruct y encoded image
fprintf('\npY...\n')
if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
    [obj.options.GAMPOptY,Op_y,EstimOut_y] = putOnGPU(obj.options.GAMPOptY,Op_y,EstimOut_y,[],[]);
end
xhatGAMP_y = gampEst_memEff(EstimIn,EstimOut_y,Op_y,obj.options.GAMPOptY);
xhatGAMP_y = clearOpts(xhatGAMP_y);
if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
    [obj.options.GAMPOptY,Op_y,EstimOut_y,~,xhatGAMP_y] = putOnCPU(obj.options.GAMPOptY,Op_y,EstimOut_y,[],xhatGAMP_y);
end

% Reconstruct z encoded image
fprintf('\npZ...\n')
if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
    [obj.options.GAMPOptZ,Op_z,EstimOut_z] = putOnGPU(obj.options.GAMPOptZ,Op_z,EstimOut_z,[],[]);
end
xhatGAMP_z = gampEst_memEff(EstimIn,EstimOut_z,Op_z,obj.options.GAMPOptZ);
xhatGAMP_z = clearOpts(xhatGAMP_z);
if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
    [obj.options.GAMPOptZ,Op_z,EstimOut_z,~,xhatGAMP_z] = putOnCPU(obj.options.GAMPOptZ,Op_z,EstimOut_z,[],xhatGAMP_z);
end

xHat_prev = [xhatGAMP_b.xhat,xhatGAMP_x.xhat,xhatGAMP_y.xhat,xhatGAMP_z.xhat];
fprintf('\nIt = %d --------------------- |dx|/|x| = %0.4f', 1,1)

% debug_save(obj,1,xhatGAMP_b,xhatGAMP_x,xhatGAMP_y,xhatGAMP_z,1,1,1)
% v_x = post_v(xhatGAMP_b.rhat, xhatGAMP_x.rhat.*maxwellCorrX, xhatGAMP_b.rvar, xhatGAMP_x.rvar, obj.options.ReVEALOpts.sigmaSq, obj.options.ReVEALOpts.gamma);
% v_y = post_v(xhatGAMP_b.rhat, xhatGAMP_y.rhat.*maxwellCorrY, xhatGAMP_b.rvar, xhatGAMP_y.rvar, obj.options.ReVEALOpts.sigmaSq, obj.options.ReVEALOpts.gamma);
% v_z = post_v(xhatGAMP_b.rhat, xhatGAMP_z.rhat.*maxwellCorrZ, xhatGAMP_b.rvar, xhatGAMP_z.rvar, obj.options.ReVEALOpts.sigmaSq, obj.options.ReVEALOpts.gamma);

%% ========================================================================
%  Estimate the Maxwell Correction from the data and apply it to the
%  ssenitvity maps
% warning('No Maxwell Correction is Applied');
% obj.data.maxwellCorrX = maxxCorr3D(ifftshift(reshape(xhatGAMP_b.xhat,size(obj.data.sampB)))...
%                         ,ifftshift(reshape(xhatGAMP_x.xhat,size(obj.data.sampB))));
% obj.data.maxwellCorrY = maxxCorr3D(ifftshift(reshape(xhatGAMP_b.xhat,size(obj.data.sampB)))...
%                         ,ifftshift(reshape(xhatGAMP_y.xhat,size(obj.data.sampB))));
% obj.data.maxwellCorrZ = maxxCorr3D(ifftshift(reshape(xhatGAMP_b.xhat,size(obj.data.sampB)))...
%                         ,ifftshift(reshape(xhatGAMP_z.xhat,size(obj.data.sampB))));


obj.data.maxwellCorrX = backgroundCorrection3D(ifftshift(ifftshift(ifftshift(reshape(xhatGAMP_b.xhat,size(obj.data.sampB)),1),2),3)...
                                              ,ifftshift(ifftshift(ifftshift(reshape(xhatGAMP_x.xhat,size(obj.data.sampB)),1),2),3));
obj.data.maxwellCorrY = backgroundCorrection3D(ifftshift(ifftshift(ifftshift(reshape(xhatGAMP_b.xhat,size(obj.data.sampB)),1),2),3)...
                                              ,ifftshift(ifftshift(ifftshift(reshape(xhatGAMP_y.xhat,size(obj.data.sampB)),1),2),3));
obj.data.maxwellCorrZ = backgroundCorrection3D(ifftshift(ifftshift(ifftshift(reshape(xhatGAMP_b.xhat,size(obj.data.sampB)),1),2),3)...
                                              ,ifftshift(ifftshift(ifftshift(reshape(xhatGAMP_z.xhat,size(obj.data.sampB)),1),2),3));

obj.data.maxwellCorrX = fftshift(fftshift(fftshift((exp(1j*obj.data.maxwellCorrX)),1),2),3);
obj.data.maxwellCorrY = fftshift(fftshift(fftshift((exp(1j*obj.data.maxwellCorrY)),1),2),3);
obj.data.maxwellCorrZ = fftshift(fftshift(fftshift((exp(1j*obj.data.maxwellCorrZ)),1),2),3);

% % apply fftshifts and restimate sensitivity maps
applyIfftshift(obj);
% I guess write out the code here b/c obj.estimateSensMaps() is not written
% in such a way to do this properly (integrate later)
%(1) Images after 1st iteration
im_b = reshape(xhatGAMP_b.xhat,size(obj.data.sampB));
im_x = reshape(xhatGAMP_x.xhat,size(obj.data.sampX));
im_y = reshape(xhatGAMP_y.xhat,size(obj.data.sampY));
im_z = reshape(xhatGAMP_z.xhat,size(obj.data.sampZ));
% (2) Apply background phase correction
im_b = im_b;
im_x = bsxfun(@times, im_x, conj(obj.data.maxwellCorrX));
im_y = bsxfun(@times, im_y, conj(obj.data.maxwellCorrY));
im_z = bsxfun(@times, im_z, conj(obj.data.maxwellCorrZ));
% (3) Multiply by current sensitivity maps
im_b = bsxfun(@times, im_b, permute(obj.data.sensMaps, [1,2,3,5,4]));
im_x = bsxfun(@times, im_x, permute(obj.data.sensMaps, [1,2,3,5,4]));
im_y = bsxfun(@times, im_y, permute(obj.data.sensMaps, [1,2,3,5,4]));
im_z = bsxfun(@times, im_z, permute(obj.data.sensMaps, [1,2,3,5,4]));
% (4) Take to k-space
k_b = fft3_shift(im_b);
k_x = fft3_shift(im_x);
k_y = fft3_shift(im_y);
k_z = fft3_shift(im_z);
% (5) Concatentate and average the kspace
k = cat(4, k_b, k_x, k_y, k_z);
k = squeeze(mean(k, 4));
% (6) Crop kspace (trust center points more)
N = size(k); M = round(size(k)/4);
k = k((N(1)/2+1)-M(1):(N(1)/2+1)+M(1), (N(2)/2+1)-M(2):(N(2)/2+1)+M(2), (N(3)/2+1)-M(3):(N(3)/2+1)+M(3), :); 

% % k(:, [1:floor(size(k,2)/4),ceil(3*(size(k,2)/4)):end], [1:floor(size(k,3)/4),ceil(3*(size(k,3)/4)):end]) = 0;
% k(:, [1:floor(size(k,2)/4),ceil(3*(size(k,2)/4)):end], :, :) = 0;
% k(:, :, [1:floor(size(k,3)/4),ceil(3*(size(k,3)/4)):end], :) = 0;
% (7) Apply windowing in kspace
% ham2(1,:,1) = hamming(ceil(3*(size(k,2)/4))-floor(size(k,2)/4));
% ham3(1,1,:) = hamming(ceil(3*(size(k,3)/4))-floor(size(k,3)/4));
ham1(:,1,1) = hamming(size(k, 1));
ham2(1,:,1) = hamming(size(k, 2));
ham3(1,1,:) = hamming(size(k, 3));
window = bsxfun(@times, ham1, bsxfun(@times, ham2, ham3));

% window = padarray(window, [0, floor(size(k,2)/4), floor(size(k,3)/4)], 'both');
k = bsxfun(@times, k, window);
M = size(k);
k = padarray(k, ceil((N-M)/2),  'pre');
k = padarray(k, floor((N-M)/2), 'post');
% (7) Bring back to image domain
im = ifft3_shift(k);
% (7) restimate sensitivity maps
p.mthd = 3; % '1' for espirit, '2' time-average spirit, '3' for walsh
p.reEst = 0; % Res-estimating sensitivities
p.fil = 3;
p.opt = 2;
[obj.data.sensMaps,~] = WalshCoilCombine3D(im,p);
obj.data.sensMaps = obj.data.sensMaps*sqrt(R);

% obj.estimateSensMaps(1);
applyFftshift(obj);

% Re-estimate sigmaSq from initialization
avg_image_b = reshape(xhatGAMP_b.xhat,size(obj.data.sampB));
avg_image_x = reshape(xhatGAMP_x.xhat,size(obj.data.sampX));
avg_image_y = reshape(xhatGAMP_y.xhat,size(obj.data.sampY));
avg_image_z = reshape(xhatGAMP_z.xhat,size(obj.data.sampZ));

avg_image = cat(4,avg_image_b,avg_image_x,avg_image_y,avg_image_z);
avg_image = mean(avg_image, 4);
avg_image = avg_image(:);
avg_image(avg_image<0.01*max(avg_image)) = [];
sigma = median(abs(avg_image));
tau = 0.15; % normally 0.15
sigmaSq = (sigma*tau)^2;
obj.options.ReVEALOpts.sigmaSq = sigmaSq*R;

obj.options.ReVEALOpts.sigmaSq = 4*obj.options.ReVEALOpts.wvar;

% Create new Input Estimation Class with Maxwell Corrections
pMRI_b = pMRI_Op_3D_t(obj.data.sensMaps,obj.data.sampB,'uniform_var',...
         obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute...
         ,'precision',obj.options.ReVEALOpts.precision,'weights',obj.data.weightsB);
pMRI_x = pMRI_Op_3D_t(bsxfun(@times,obj.data.sensMaps,obj.data.maxwellCorrX),obj.data.sampX,'uniform_var',...
         obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute...
         ,'precision',obj.options.ReVEALOpts.precision,'weights',obj.data.weightsX);
pMRI_y = pMRI_Op_3D_t(bsxfun(@times,obj.data.sensMaps,obj.data.maxwellCorrY),obj.data.sampY,'uniform_var',...
         obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute...
         ,'precision',obj.options.ReVEALOpts.precision,'weights',obj.data.weightsY);
pMRI_z = pMRI_Op_3D_t(bsxfun(@times,obj.data.sensMaps,obj.data.maxwellCorrZ),obj.data.sampZ,'uniform_var',...
         obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute...
         ,'precision',obj.options.ReVEALOpts.precision,'weights',obj.data.weightsZ);

% Create nd-DWT Linear Transform
W = ndDWTLinTrans4D(obj.options.SparseTrans,size(obj.data.sampB),'uniform_var',...
    obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute,...
    'precision',obj.options.ReVEALOpts.precision);

% Concatonate All The Linear Transform Operator Together
Op_b = LinTransConcat({pMRI_b;W},[1,1],obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
Op_x = LinTransConcat({pMRI_x;W},[1,1],obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
Op_y = LinTransConcat({pMRI_y;W},[1,1],obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
Op_z = LinTransConcat({pMRI_z;W},[1,1],obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);

clear pMRI_b pMRI_x pMRI_y pMRI_z W dataB dataX dataY dataZ

debug_save(obj,1,xhatGAMP_b,xhatGAMP_x,xhatGAMP_y,xhatGAMP_z,ones(size(xhatGAMP_b.xhat)),ones(size(xhatGAMP_b.xhat)),ones(size(xhatGAMP_b.xhat)))

%% =========================================================================
% First-order Markov chain properties
p01 = 0.05;     % Default active-to-inactive transition probability
% learn_p01 = 'true';      % Learn p01 using EM alg. by default
dim = 'row';    % Each row forms a Markov chain
sparsity = 0.2;
N = size(obj.data.sampB);
sparsity = sparsity*ones([N(1)*N(2)*N(3),N(4)]);
mrc = MarkovChain1('p01',p01,'dim',dim);

% tmp = reshape(obj.data.x0,size(obj.data.sampB));
% level = graythresh(abs(tmp(:,:,1)));
% mask = im2bw(abs(tmp(:,:,1)),level);
% mask = repmat(mask,[1,1,size(obj.data.sampB,3)]);
% clear tmp

% % Calculate posterior on velocity support
% gamma = obj.options.ReVEALOpts.gamma;
% v_x = post_v(xhatGAMP_b.rhat, xhatGAMP_x.rhat, xhatGAMP_b.rvar, xhatGAMP_x.rvar, obj.options.ReVEALOpts.sigmaSq, gamma);
% v_y = post_v(xhatGAMP_b.rhat, xhatGAMP_y.rhat, xhatGAMP_b.rvar, xhatGAMP_y.rvar, obj.options.ReVEALOpts.sigmaSq, gamma);
% v_z = post_v(xhatGAMP_b.rhat, xhatGAMP_z.rhat, xhatGAMP_b.rvar, xhatGAMP_z.rvar, obj.options.ReVEALOpts.sigmaSq, gamma);
% 
% % Take the posterior over all channels
% v = max([v_x(:),v_y(:),v_z(:)],[],2);
% v(v<1e-6) = 1e-6;
% v(v>=1) = 0.999;
% N = size(obj.data.sampB);
% v = reshape(v,[N(1)*N(2)*N(3),N(4)]);
% 
% % Take the max or MEAN? over time
% % obj.options.ReVEALOpts.gamma = max(reshape(v,size(obj.data.sampB)),[],4);
% obj.options.ReVEALOpts.gamma = mean(reshape(v, size(obj.data.sampB)), 4);
% obj.options.ReVEALOpts.gamma = repmat(obj.options.ReVEALOpts.gamma,[1,1,1,size(obj.data.sampB,4)]);
% obj.options.ReVEALOpts.gamma = obj.options.ReVEALOpts.gamma(:);
% gamma = obj.options.ReVEALOpts.gamma;

gamma = obj.options.ReVEALOpts.gamma;

%% =========================================================================
% Main Loop
for ind = 2:obj.options.ReVEALOpts.nit

    obj.options.GAMPOptB.nit = 10;
    obj.options.GAMPOptX.nit = 10;
    obj.options.GAMPOptY.nit = 10;
    obj.options.GAMPOptZ.nit = 10;
       
    %======================================================================
    % Create new input estimators
    pvarb = xhatGAMP_b.rvar + obj.options.ReVEALOpts.sigmaSq;
    pvarx = xhatGAMP_x.rvar + obj.options.ReVEALOpts.sigmaSq;
    pvary = xhatGAMP_y.rvar + obj.options.ReVEALOpts.sigmaSq;
    pvarz = xhatGAMP_z.rvar + obj.options.ReVEALOpts.sigmaSq;
    
    % Input Estimator for pb
    meanB = (abs(xhatGAMP_x.rhat).*pvary.*pvarz +...
             abs(xhatGAMP_y.rhat).*pvarx.*pvarz +...
             abs(xhatGAMP_z.rhat).*pvarx.*pvary)./...
             (pvarx.*pvary + pvarx.*pvarz + pvary.*pvarz) + eps;
    varB = (pvarx.*pvary.*pvarz)./...
             (pvarx.*pvary + pvarx.*pvarz + pvary.*pvarz) + eps;
    EstimIn_b1 = ncCAwgnEstimIn(meanB,varB);
    meanB = (xhatGAMP_x.rhat.*pvary.*pvarz +...
             xhatGAMP_y.rhat.*pvarx.*pvarz +...
             xhatGAMP_z.rhat.*pvarx.*pvary)./...
             (pvarx.*pvary + pvarx.*pvarz + pvary.*pvarz) + eps;
    EstimIn_b2 = CAwgnEstimIn(meanB,varB);
    EstimIn_b = MixScaEstimIn(EstimIn_b1,obj.options.ReVEALOpts.gamma,EstimIn_b2);

    % Input Estimator for px
    meanX = (abs(xhatGAMP_b.rhat).*pvary.*pvarz +...
             abs(xhatGAMP_y.rhat).*pvarb.*pvarz +...
             abs(xhatGAMP_z.rhat).*pvarb.*pvary)./...
             (pvarb.*pvary + pvarb.*pvarz + pvary.*pvarz) + eps;
    varX = (pvarb.*pvary.*pvarz)./...
             (pvarb.*pvary + pvarb.*pvarz + pvary.*pvarz) + eps;     
    EstimIn_x1 = ncCAwgnEstimIn(meanX,varX);
    meanX = (xhatGAMP_b.rhat.*pvary.*pvarz +...
             xhatGAMP_y.rhat.*pvarb.*pvarz +...
             xhatGAMP_z.rhat.*pvarb.*pvary)./...
             (pvarb.*pvary + pvarb.*pvarz + pvary.*pvarz) + eps;
    EstimIn_x2 = CAwgnEstimIn(meanX,varX);
    EstimIn_x = MixScaEstimIn(EstimIn_x1,obj.options.ReVEALOpts.gamma,EstimIn_x2);

    % Input Estimator for py
    meanY = (abs(xhatGAMP_b.rhat).*pvarx.*pvarz +...
             abs(xhatGAMP_x.rhat).*pvarb.*pvarz +...
             abs(xhatGAMP_z.rhat).*pvarb.*pvarx)./...
             (pvarb.*pvarx + pvarb.*pvarz + pvarx.*pvarz) + eps;
    varY = (pvarb.*pvarx.*pvarz)./...
             (pvarb.*pvarx + pvarb.*pvarz + pvarx.*pvarz) + eps;     
    EstimIn_y1 = ncCAwgnEstimIn(meanY,varY);
    meanY = (xhatGAMP_b.rhat.*pvarx.*pvarz +...
             xhatGAMP_x.rhat.*pvarb.*pvarz +...
             xhatGAMP_z.rhat.*pvarb.*pvarx)./...
             (pvarb.*pvarx + pvarb.*pvarz + pvarx.*pvarz) + eps;
    EstimIn_y2 = CAwgnEstimIn(meanY,varY);
    EstimIn_y = MixScaEstimIn(EstimIn_y1,obj.options.ReVEALOpts.gamma,EstimIn_y2);

    % Input Estimator for pz
    meanZ = (abs(xhatGAMP_b.rhat).*pvarx.*pvary +...
             abs(xhatGAMP_x.rhat).*pvarb.*pvary +...
             abs(xhatGAMP_y.rhat).*pvarb.*pvarx)./...
             (pvarb.*pvarx + pvarb.*pvary + pvarx.*pvary) + eps;
    varZ = (pvarb.*pvarx.*pvary)./...
             (pvarb.*pvarx + pvarb.*pvary + pvarx.*pvary) + eps;
    EstimIn_z1 = ncCAwgnEstimIn(meanZ,varZ);
    meanZ = (xhatGAMP_b.rhat.*pvarx.*pvary +...
             xhatGAMP_x.rhat.*pvarb.*pvary +...
             xhatGAMP_y.rhat.*pvarb.*pvarx)./...
             (pvarb.*pvarx + pvarb.*pvary + pvarx.*pvary) + eps;
    EstimIn_z2 = CAwgnEstimIn(meanZ,varZ);
    EstimIn_z = MixScaEstimIn(EstimIn_z1,obj.options.ReVEALOpts.gamma,EstimIn_z2);
    
    % Clear some memory up
    clear meanB meanX meanY meanZ varB varX varY varZ pvarb pvarx pvary...
          pvarz EstimIn_b1 EstimIn_b2 EstimIn_x1 EstimIn_x2 EstimIn_y1...
          EstimIn_y2 EstimIn_z1 EstimIn_z2;

    %======================================================================
    % Warm Start xb
    obj.options.GAMPOptB.xhat0 = xhatGAMP_b.xhat;
    obj.options.GAMPOptB.xvar0 = xhatGAMP_b.xvar;
    obj.options.GAMPOptB.shat0 = xhatGAMP_b.shat;
    obj.options.GAMPOptB.svar0 = xhatGAMP_b.svar;
    obj.options.GAMPOptB.xhatPrev0 = xhatGAMP_b.xhatPrev;
    obj.options.GAMPOptB.scaleFac = xhatGAMP_b.scaleFac;
    obj.options.GAMPOptB.step = min(max(xhatGAMP_b.step,0.05),xhatGAMP_b.stepMax);

    % Warm Start xx
    obj.options.GAMPOptX.xhat0 = xhatGAMP_x.xhat;
    obj.options.GAMPOptX.xvar0 = xhatGAMP_x.xvar;
    obj.options.GAMPOptX.shat0 = xhatGAMP_x.shat;
    obj.options.GAMPOptX.svar0 = xhatGAMP_x.svar;
    obj.options.GAMPOptX.xhatPrev0 = xhatGAMP_x.xhatPrev;
    obj.options.GAMPOptX.scaleFac = xhatGAMP_x.scaleFac;
    obj.options.GAMPOptX.step = min(max(xhatGAMP_x.step,0.05),xhatGAMP_x.stepMax);

    % Warm Start xy
    obj.options.GAMPOptY.xhat0 = xhatGAMP_y.xhat;
    obj.options.GAMPOptY.xvar0 = xhatGAMP_y.xvar;
    obj.options.GAMPOptY.shat0 = xhatGAMP_y.shat;
    obj.options.GAMPOptY.svar0 = xhatGAMP_y.svar;
    obj.options.GAMPOptY.xhatPrev0 = xhatGAMP_y.xhatPrev;
    obj.options.GAMPOptY.scaleFac = xhatGAMP_y.scaleFac;
    obj.options.GAMPOptY.step = min(max(xhatGAMP_y.step,0.05),xhatGAMP_y.stepMax);

    % Warm Start xz
    obj.options.GAMPOptZ.xhat0 = xhatGAMP_z.xhat;
    obj.options.GAMPOptZ.xvar0 = xhatGAMP_z.xvar;
    obj.options.GAMPOptZ.shat0 = xhatGAMP_z.shat;
    obj.options.GAMPOptZ.svar0 = xhatGAMP_z.svar;
    obj.options.GAMPOptZ.xhatPrev0 = xhatGAMP_z.xhatPrev;
    obj.options.GAMPOptZ.scaleFac = xhatGAMP_z.scaleFac;
    obj.options.GAMPOptZ.step = min(max(xhatGAMP_z.step,0.05),xhatGAMP_z.stepMax);
    
    % =====================================================================
    % MRC
    if ind >40
        % Calculate Posterior on velocity support
        v_x = post_v(xhatGAMP_b.rhat, xhatGAMP_x.rhat, xhatGAMP_b.rvar, xhatGAMP_x.rvar, obj.options.ReVEALOpts.sigmaSq, gamma);
        v_y = post_v(xhatGAMP_b.rhat, xhatGAMP_y.rhat, xhatGAMP_b.rvar, xhatGAMP_y.rvar, obj.options.ReVEALOpts.sigmaSq, gamma);
        v_z = post_v(xhatGAMP_b.rhat, xhatGAMP_z.rhat, xhatGAMP_b.rvar, xhatGAMP_z.rvar, obj.options.ReVEALOpts.sigmaSq, gamma);
        
        % Take the posterior over all channels
        v = max([v_x(:),v_y(:),v_z(:)],[],2);
        v(v<1e-6) = 1e-6;
        v(v>=1) = 0.999;
        N = size(obj.data.sampB);
        v = reshape(v,[N(1)*N(2)*N(3),N(4)]);
        
        % Update Chain
        [v, S_POST, p01_upd] = mrc.binarymarkov_msgs(v, sparsity);
        %sparsity = mean(S_POST(mask ==1));
        sparsity = mean(S_POST(:));
        sparsity = sparsity*ones([N(1)*N(2)*N(3),N(4)]);
        mrc.p01 = p01_upd;
        
        % Take the max over time
        obj.options.ReVEALOpts.gamma = max(reshape(v,size(obj.data.sampB)),[],4);
        obj.options.ReVEALOpts.gamma = repmat(obj.options.ReVEALOpts.gamma,[1,1,1,size(obj.data.sampB,4)]);
        obj.options.ReVEALOpts.gamma = obj.options.ReVEALOpts.gamma(:);
        gamma = v(:);
        fprintf(sprintf('\nEM Update p01 = %s, sparsity = %s\n',num2str(p01_upd),num2str(sparsity(1))));
    end
    
%     % Calculate posterior on velocity support
%     v_x = post_v(xhatGAMP_b.rhat, xhatGAMP_x.rhat, xhatGAMP_b.rvar, xhatGAMP_x.rvar, obj.options.ReVEALOpts.sigmaSq, gamma);
%     v_y = post_v(xhatGAMP_b.rhat, xhatGAMP_y.rhat, xhatGAMP_b.rvar, xhatGAMP_y.rvar, obj.options.ReVEALOpts.sigmaSq, gamma);
%     v_z = post_v(xhatGAMP_b.rhat, xhatGAMP_z.rhat, xhatGAMP_b.rvar, xhatGAMP_z.rvar, obj.options.ReVEALOpts.sigmaSq, gamma);
%     
%     % Take the posterior over all channels
%     v = max([v_x(:),v_y(:),v_z(:)],[],2); % max() here might be a problem
%     v(v<1e-6) = 1e-6;
%     v(v>=1) = 0.999;
%     N = size(obj.data.sampB);
%     v = reshape(v,[N(1)*N(2)*N(3),N(4)]);
%     
%     % Take the max over time
% %     obj.options.ReVEALOpts.gamma = max(reshape(v,size(obj.data.sampB)),[],4);
%     obj.options.ReVEALOpts.gamma = mean(reshape(v, size(obj.data.sampB)), 4);
%     obj.options.ReVEALOpts.gamma = repmat(obj.options.ReVEALOpts.gamma,[1,1,1,size(obj.data.sampB,4)]);
%     obj.options.ReVEALOpts.gamma = obj.options.ReVEALOpts.gamma(:);
%     gamma = obj.options.ReVEALOpts.gamma;
    
    % Clear some memory up
    clear v xhatGAMP_b xhatGAMP_x xhatGAMP_y xhatGAMP_z; 

    %======================================================================
    % Reconstruct
    fprintf('\npB...\n')
    if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
        [obj.options.GAMPOptB,Op_b,EstimOut_b] = putOnGPU(obj.options.GAMPOptB,Op_b,EstimOut_b,EstimIn_b,[]);
    end
    xhatGAMP_b = gampEst_memEff(EstimIn_b,EstimOut_b,Op_b,obj.options.GAMPOptB);
    xhatGAMP_b = clearOpts(xhatGAMP_b);
    if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
        [obj.options.GAMPOptB,Op_b,EstimOut_b,~,xhatGAMP_b] = putOnCPU(obj.options.GAMPOptB,Op_b,EstimOut_b,EstimIn_b,xhatGAMP_b);
    end

    % Reconstruct x encoded image
    fprintf('\npX...\n')
    if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
        [obj.options.GAMPOptX,Op_x,EstimOut_x] = putOnGPU(obj.options.GAMPOptX,Op_x,EstimOut_x,EstimIn_x,[]);
    end
    xhatGAMP_x = gampEst_memEff(EstimIn_x,EstimOut_x,Op_x,obj.options.GAMPOptX);
    xhatGAMP_x = clearOpts(xhatGAMP_x);
    if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
        [obj.options.GAMPOptX,Op_x,EstimOut_x,~,xhatGAMP_x] = putOnCPU(obj.options.GAMPOptX,Op_x,EstimOut_x,EstimIn_x,xhatGAMP_x);
    end

    % Reconstruct y encoded image
    fprintf('\npY...\n')
    if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
        [obj.options.GAMPOptY,Op_y,EstimOut_y] = putOnGPU(obj.options.GAMPOptY,Op_y,EstimOut_y,EstimIn_y,[]);
    end
    xhatGAMP_y = gampEst_memEff(EstimIn_y,EstimOut_y,Op_y,obj.options.GAMPOptY);
    xhatGAMP_y = clearOpts(xhatGAMP_y);
    if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
        [obj.options.GAMPOptY,Op_y,EstimOut_y,~,xhatGAMP_y] = putOnCPU(obj.options.GAMPOptY,Op_y,EstimOut_y,EstimIn_y,xhatGAMP_y);
    end

    % Reconstruct z encoded image
    fprintf('\npZ...\n')
    if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
        [obj.options.GAMPOptZ,Op_z,EstimOut_z] = putOnGPU(obj.options.GAMPOptZ,Op_z,EstimOut_z,EstimIn_z,[]);
    end
    xhatGAMP_z = gampEst_memEff(EstimIn_z,EstimOut_z,Op_z,obj.options.GAMPOptZ);
    xhatGAMP_z = clearOpts(xhatGAMP_z);
    if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
        [obj.options.GAMPOptZ,Op_z,EstimOut_z,~,xhatGAMP_z] = putOnCPU(obj.options.GAMPOptZ,Op_z,EstimOut_z,EstimIn_z,xhatGAMP_z);
    end
    
    
    x = [xhatGAMP_b.xhat,xhatGAMP_x.xhat,xhatGAMP_y.xhat,xhatGAMP_z.xhat];
    dx= norm(xHat_prev-x)/norm(x);
    xHat_prev = x;
    clear x;
    fprintf('\nIt = %d --------------------- |dx|/|x| = %0.4f', ind,dx)
    debug_save(obj,ind+1,xhatGAMP_b,xhatGAMP_x,xhatGAMP_y,xhatGAMP_z,ones(size(xhatGAMP_b.xhat)),ones(size(xhatGAMP_b.xhat)),ones(size(xhatGAMP_b.xhat)))
    
%     if dx < 1e-2 % less than 1% change
%         disp('Terminating iterations early:');
%         disp(['Less than ',num2str(100*1e-2),'% change in x between iterations.']);
%         break;
%     end
%     debug_save(obj,ind,xhatGAMP_b,xhatGAMP_x,xhatGAMP_y,xhatGAMP_z,obj.data.maxwellCorrX,obj.data.maxwellCorrY,obj.data.maxwellCorrZ)
end
fprintf('ReVAMP Calculated in %s s\n', num2str(toc))

%% Calculate velocity maps
obj.outputs.thetaX = angle(xhatGAMP_x.xhat.*conj(xhatGAMP_b.xhat));
obj.outputs.thetaY = angle(xhatGAMP_y.xhat.*conj(xhatGAMP_b.xhat));
obj.outputs.thetaZ = angle(xhatGAMP_z.xhat.*conj(xhatGAMP_b.xhat));

v_x = post_v(xhatGAMP_b.rhat, xhatGAMP_x.rhat, xhatGAMP_b.rvar, xhatGAMP_x.rvar, obj.options.ReVEALOpts.sigmaSq, obj.options.ReVEALOpts.gamma);
v_y = post_v(xhatGAMP_b.rhat, xhatGAMP_y.rhat, xhatGAMP_b.rvar, xhatGAMP_y.rvar, obj.options.ReVEALOpts.sigmaSq, obj.options.ReVEALOpts.gamma);
v_z = post_v(xhatGAMP_b.rhat, xhatGAMP_z.rhat, xhatGAMP_b.rvar, xhatGAMP_z.rvar, obj.options.ReVEALOpts.sigmaSq, obj.options.ReVEALOpts.gamma);

% Reshape outputs
obj.outputs.xHatb = reshape(xhatGAMP_b.xhat,size(obj.data.sampB));
obj.outputs.xHatx = reshape(xhatGAMP_x.xhat,size(obj.data.sampX));
obj.outputs.xHaty = reshape(xhatGAMP_y.xhat,size(obj.data.sampY));
obj.outputs.xHatz = reshape(xhatGAMP_z.xhat,size(obj.data.sampZ));

obj.outputs.thetaX = reshape(obj.outputs.thetaX,size(obj.data.sampB));
obj.outputs.thetaY = reshape(obj.outputs.thetaY,size(obj.data.sampB));
obj.outputs.thetaZ = reshape(obj.outputs.thetaZ,size(obj.data.sampB));

obj.outputs.vX = reshape(v_x,size(obj.data.sampB));
obj.outputs.vY = reshape(v_y,size(obj.data.sampB));
obj.outputs.vZ = reshape(v_z,size(obj.data.sampB));

obj.outputs.scanParam = obj.data.scanParam;

%% fftshifts
% warning('Change to ifftshift')
obj.data.Yb = ifftshift(ifftshift(ifftshift(obj.data.Yb,1),2),3);
obj.data.Yx = ifftshift(ifftshift(ifftshift(obj.data.Yx,1),2),3);
obj.data.Yy = ifftshift(ifftshift(ifftshift(obj.data.Yy,1),2),3);
obj.data.Yz = ifftshift(ifftshift(ifftshift(obj.data.Yz,1),2),3);

obj.data.sampB = ifftshift(ifftshift(ifftshift(obj.data.sampB,1),2),3);
obj.data.sampX = ifftshift(ifftshift(ifftshift(obj.data.sampX,1),2),3);
obj.data.sampY = ifftshift(ifftshift(ifftshift(obj.data.sampY,1),2),3);
obj.data.sampZ = ifftshift(ifftshift(ifftshift(obj.data.sampZ,1),2),3);

obj.data.maxwellCorrX =  ifftshift(ifftshift(ifftshift(obj.data.maxwellCorrX,1),2),3);
obj.data.maxwellCorrY =  ifftshift(ifftshift(ifftshift(obj.data.maxwellCorrY,1),2),3);
obj.data.maxwellCorrZ =  ifftshift(ifftshift(ifftshift(obj.data.maxwellCorrZ,1),2),3);

obj.data.sensMaps = fftshift(fftshift(fftshift(obj.data.sensMaps,1),2),3);
obj.data.x0 =  fftshift(fftshift(fftshift(obj.data.x0,1),2),3);
% obj.options.ReVEALOpts.gamma = fftshift(fftshift(fftshift(obj.options.ReVEALOpts.gamma,1),2),3);
% obj.options.ReVEALOpts.gamma = obj.options.ReVEALOpts.gamma(:);

obj.outputs.xHatb = ifftshift(ifftshift(ifftshift(obj.outputs.xHatb,1),2),3);
obj.outputs.xHatx = ifftshift(ifftshift(ifftshift(obj.outputs.xHatx,1),2),3);
obj.outputs.xHaty = ifftshift(ifftshift(ifftshift(obj.outputs.xHaty,1),2),3);
obj.outputs.xHatz = ifftshift(ifftshift(ifftshift(obj.outputs.xHatz,1),2),3);

obj.outputs.thetaX = ifftshift(ifftshift(ifftshift(obj.outputs.thetaX,1),2),3);
obj.outputs.thetaY = ifftshift(ifftshift(ifftshift(obj.outputs.thetaY,1),2),3);
obj.outputs.thetaZ = ifftshift(ifftshift(ifftshift(obj.outputs.thetaZ,1),2),3);

obj.outputs.vX = ifftshift(ifftshift(ifftshift(obj.outputs.vX,1),2),3);
obj.outputs.vY = ifftshift(ifftshift(ifftshift(obj.outputs.vY,1),2),3);
obj.outputs.vZ = ifftshift(ifftshift(ifftshift(obj.outputs.vZ,1),2),3);


end

% Put data on the GPU
function [GAMPOpt,Op,outputEstim,inputEstim,GAMPOutput] = putOnGPU(GAMPOpt,Op,outputEstim,inputEstim,GAMPOutput)
    GAMPOpt.xhat0 = gpuArray(GAMPOpt.xhat0);

    if ~isempty(GAMPOutput)
        GAMPOutput.xhat = gpuArray(GAMPOutput.xhat);
        GAMPOutput.xvar = gpuArray(GAMPOutput.xvar);
        GAMPOutput.rhat = gpuArray(GAMPOutput.rhat);
        GAMPOutput.rvar = gpuArray(GAMPOutput.rvar);
    end
    
    Op.lta{1}.C = gpuArray(Op.lta{1}.C);
    Op.lta{1}.CSq= gpuArray(Op.lta{1}.CSq);
    Op.lta{1}.CSqTr = gpuArray(Op.lta{1}.CSqTr);
    Op.lta{1}.uniform_var = gpuArray(Op.lta{1}.uniform_var);
    Op.lta{1}.mask_patterns = gpuArray(Op.lta{1}.mask_patterns);
    Op.lta{1}.mask_weights = gpuArray(Op.lta{1}.mask_weights);
    
    
    outputEstim.estimArray{1}.y = gpuArray(outputEstim.estimArray{1}.y);
    
    if ~isempty(inputEstim)
        inputEstim.estim1.mean0 = gpuArray(inputEstim.estim1.mean0);
        inputEstim.estim1.var0 = gpuArray(inputEstim.estim1.var0);
        inputEstim.estim0.mean0 = gpuArray(inputEstim.estim0.mean0);
        inputEstim.estim0.var0 = gpuArray(inputEstim.estim0.var0);
    end
end

% Put data on the CPU
function [GAMPOpt,Op,outputEstim,inputEstim,GAMPOutput] = putOnCPU(GAMPOpt,Op,outputEstim,inputEstim,GAMPOutput)
    GAMPOpt.xhat0 = gather(GAMPOpt.xhat0);

    if ~isempty(GAMPOutput)
        GAMPOutput.xhat = gather(GAMPOutput.xhat);
        GAMPOutput.xvar = gather(GAMPOutput.xvar);
        GAMPOutput.rhat = gather(GAMPOutput.rhat);
        GAMPOutput.rvar = gather(GAMPOutput.rvar);
    end
    
    Op.lta{1}.C = gather(Op.lta{1}.C);
    Op.lta{1}.CSq= gather(Op.lta{1}.CSq);
    Op.lta{1}.CSqTr = gather(Op.lta{1}.CSqTr);
    Op.lta{1}.uniform_var = gather(Op.lta{1}.uniform_var);
    Op.lta{1}.mask_patterns = gather(Op.lta{1}.mask_patterns);
    Op.lta{1}.mask_weights = gather(Op.lta{1}.mask_weights);
    
    outputEstim.estimArray{1}.y = gather(outputEstim.estimArray{1}.y);
    
    if ~isempty(inputEstim)
        inputEstim.estim1.mean0 = gather(inputEstim.estim1.mean0);
        inputEstim.estim1.var0 = gather(inputEstim.estim1.var0);
        inputEstim.estim0.mean0 = gather(inputEstim.estim0.mean0);
        inputEstim.estim0.var0 = gather(inputEstim.estim0.var0);
    end
end

% Clear out unneeded GAMP outputs to save memory
function GAMPOpt = clearOpts(GAMPOpt)
    %GAMPOpt.xhat
    %GAMPOpt.xvar
    GAMPOpt.phat = [];
    GAMPOpt.pvar = [];
    GAMPOpt.zhat = [];
    GAMPOpt.zvar = [];
    GAMPOpt.shat = [];
    GAMPOpt.svar = [];
    %GAMPOpt.rhat
    %GAMPOpt.rvar
    GAMPOpt.Axhat = [];
    GAMPOpt.xhatPrev = [];
    GAMPOpt.xhatNext = [];
    GAMPOpt.xvarNext = [];
    GAMPOpt.xhatDamp = [];
    GAMPOpt.pvarOpt = [];
    GAMPOpt.rvarOpt = [];
    GAMPOpt.A2xvarOpt = [];
    GAMPOpt.shatNext = [];
    GAMPOpt.svarNext = [];
end

function debug_save(obj,it,xhatGAMP_b,xhatGAMP_x,xhatGAMP_y,xhatGAMP_z,maxwellCorrX,maxwellCorrY,maxwellCorrZ)
    %% Calculate velocity maps
    maxwellCorrX = 1;
    maxwellCorrY = 1;
    maxwellCorrZ = 1;
    obj.outputs.thetaX = angle(xhatGAMP_x.xhat.*conj(xhatGAMP_b.xhat).*maxwellCorrX);
    obj.outputs.thetaY = angle(xhatGAMP_y.xhat.*conj(xhatGAMP_b.xhat).*maxwellCorrY);
    obj.outputs.thetaZ = angle(xhatGAMP_z.xhat.*conj(xhatGAMP_b.xhat).*maxwellCorrZ);

    v_x = post_v(xhatGAMP_b.rhat, xhatGAMP_x.rhat.*maxwellCorrX, xhatGAMP_b.rvar, xhatGAMP_x.rvar, obj.options.ReVEALOpts.sigmaSq, obj.options.ReVEALOpts.gamma);
    v_y = post_v(xhatGAMP_b.rhat, xhatGAMP_y.rhat.*maxwellCorrY, xhatGAMP_b.rvar, xhatGAMP_y.rvar, obj.options.ReVEALOpts.sigmaSq, obj.options.ReVEALOpts.gamma);
    v_z = post_v(xhatGAMP_b.rhat, xhatGAMP_z.rhat.*maxwellCorrZ, xhatGAMP_b.rvar, xhatGAMP_z.rvar, obj.options.ReVEALOpts.sigmaSq, obj.options.ReVEALOpts.gamma);

    % Reshape outputs
    obj.outputs.xHatb = reshape(xhatGAMP_b.xhat,size(obj.data.sampB));
    obj.outputs.xHatx = reshape(xhatGAMP_x.xhat,size(obj.data.sampX));
    obj.outputs.xHaty = reshape(xhatGAMP_y.xhat,size(obj.data.sampY));
    obj.outputs.xHatz = reshape(xhatGAMP_z.xhat,size(obj.data.sampZ));

    obj.outputs.thetaX = reshape(obj.outputs.thetaX,size(obj.data.sampB));
    obj.outputs.thetaY = reshape(obj.outputs.thetaY,size(obj.data.sampB));
    obj.outputs.thetaZ = reshape(obj.outputs.thetaZ,size(obj.data.sampB));

    obj.outputs.vX = reshape(v_x,size(obj.data.sampB));
    obj.outputs.vY = reshape(v_y,size(obj.data.sampB));
    obj.outputs.vZ = reshape(v_z,size(obj.data.sampB));

    %% fftshifts
    warning('Change to ifftshift')
    obj.data.Yb = fftshift(fftshift(fftshift(obj.data.Yb,1),2),3);
    obj.data.Yx = fftshift(fftshift(fftshift(obj.data.Yx,1),2),3);
    obj.data.Yy = fftshift(fftshift(fftshift(obj.data.Yy,1),2),3);
    obj.data.Yz = fftshift(fftshift(fftshift(obj.data.Yz,1),2),3);

    obj.data.sampB = fftshift(fftshift(fftshift(obj.data.sampB,1),2),3);
    obj.data.sampX = fftshift(fftshift(fftshift(obj.data.sampX,1),2),3);
    obj.data.sampY = fftshift(fftshift(fftshift(obj.data.sampY,1),2),3);
    obj.data.sampZ = fftshift(fftshift(fftshift(obj.data.sampZ,1),2),3);

%     obj.data.maxwellCorrX =  fftshift(fftshift(fftshift(obj.data.maxwellCorrX,1),2),3);
%     obj.data.maxwellCorrY =  fftshift(fftshift(fftshift(obj.data.maxwellCorrY,1),2),3);
%     obj.data.maxwellCorrZ =  fftshift(fftshift(fftshift(obj.data.maxwellCorrZ,1),2),3);

    obj.data.sensMaps = fftshift(fftshift(fftshift(obj.data.sensMaps,1),2),3);
    obj.data.x0 =  fftshift(fftshift(fftshift(obj.data.x0,1),2),3);
    % obj.options.ReVEALOpts.gamma = fftshift(fftshift(fftshift(obj.options.ReVEALOpts.gamma,1),2),3);
    % obj.options.ReVEALOpts.gamma = obj.options.ReVEALOpts.gamma(:);

    obj.outputs.xHatb = fftshift(fftshift(fftshift(obj.outputs.xHatb,1),2),3);
    obj.outputs.xHatx = fftshift(fftshift(fftshift(obj.outputs.xHatx,1),2),3);
    obj.outputs.xHaty = fftshift(fftshift(fftshift(obj.outputs.xHaty,1),2),3);
    obj.outputs.xHatz = fftshift(fftshift(fftshift(obj.outputs.xHatz,1),2),3);

    obj.outputs.thetaX = fftshift(fftshift(fftshift(obj.outputs.thetaX,1),2),3);
    obj.outputs.thetaY = fftshift(fftshift(fftshift(obj.outputs.thetaY,1),2),3);
    obj.outputs.thetaZ = fftshift(fftshift(fftshift(obj.outputs.thetaZ,1),2),3);

    obj.outputs.vX = fftshift(fftshift(fftshift(obj.outputs.vX,1),2),3);
    obj.outputs.vY = fftshift(fftshift(fftshift(obj.outputs.vY,1),2),3);
    obj.outputs.vZ = fftshift(fftshift(fftshift(obj.outputs.vZ,1),2),3);
    
    lambda = obj.options.ReVEALOpts.lambda0;
    sigma_sq = obj.options.ReVEALOpts.sigmaSq;
    
    name = ['ReVEAL_Iter_',num2str(it),'_lambda_',num2str(lambda),'_sigma_sq_',num2str(sigma_sq)];
    
    if ~exist(name,'dir')
        mkdir(name);
    end
        
    obj.saveGIFs([name,'\',name])
end

% apply an fftshift to the input data
function applyFftshift(obj)
    obj.data.Yb = fftshift(fftshift(fftshift(obj.data.Yb,1),2),3);
    obj.data.Yx = fftshift(fftshift(fftshift(obj.data.Yx,1),2),3);
    obj.data.Yy = fftshift(fftshift(fftshift(obj.data.Yy,1),2),3);
    obj.data.Yz = fftshift(fftshift(fftshift(obj.data.Yz,1),2),3);

    obj.data.sampB = fftshift(fftshift(fftshift(obj.data.sampB,1),2),3);
    obj.data.sampX = fftshift(fftshift(fftshift(obj.data.sampX,1),2),3);
    obj.data.sampY = fftshift(fftshift(fftshift(obj.data.sampY,1),2),3);
    obj.data.sampZ = fftshift(fftshift(fftshift(obj.data.sampZ,1),2),3);
    
    obj.data.weightsB = fftshift(fftshift(fftshift(obj.data.weightsB,1),2),3);
    obj.data.weightsX = fftshift(fftshift(fftshift(obj.data.weightsX,1),2),3);
    obj.data.weightsY = fftshift(fftshift(fftshift(obj.data.weightsY,1),2),3);
    obj.data.weightsZ = fftshift(fftshift(fftshift(obj.data.weightsZ,1),2),3);
    
    obj.data.maxwellCorrX =  fftshift(fftshift(fftshift(obj.data.maxwellCorrX,1),2),3);
    obj.data.maxwellCorrY =  fftshift(fftshift(fftshift(obj.data.maxwellCorrY,1),2),3);
    obj.data.maxwellCorrZ =  fftshift(fftshift(fftshift(obj.data.maxwellCorrZ,1),2),3);

    obj.data.sensMaps = fftshift(fftshift(fftshift(obj.data.sensMaps,1),2),3);
    obj.data.x0 =  fftshift(fftshift(fftshift(obj.data.x0,1),2),3);
end

% apply an ifftshift to the input data
function applyIfftshift(obj)
    obj.data.Yb = ifftshift(ifftshift(ifftshift(obj.data.Yb,1),2),3);
    obj.data.Yx = ifftshift(ifftshift(ifftshift(obj.data.Yx,1),2),3);
    obj.data.Yy = ifftshift(ifftshift(ifftshift(obj.data.Yy,1),2),3);
    obj.data.Yz = ifftshift(ifftshift(ifftshift(obj.data.Yz,1),2),3);

    obj.data.sampB = ifftshift(ifftshift(ifftshift(obj.data.sampB,1),2),3);
    obj.data.sampX = ifftshift(ifftshift(ifftshift(obj.data.sampX,1),2),3);
    obj.data.sampY = ifftshift(ifftshift(ifftshift(obj.data.sampY,1),2),3);
    obj.data.sampZ = ifftshift(ifftshift(ifftshift(obj.data.sampZ,1),2),3);
    
    obj.data.weightsB = ifftshift(ifftshift(ifftshift(obj.data.weightsB,1),2),3);
    obj.data.weightsX = ifftshift(ifftshift(ifftshift(obj.data.weightsX,1),2),3);
    obj.data.weightsY = ifftshift(ifftshift(ifftshift(obj.data.weightsY,1),2),3);
    obj.data.weightsZ = ifftshift(ifftshift(ifftshift(obj.data.weightsZ,1),2),3);

    obj.data.maxwellCorrX =  ifftshift(ifftshift(ifftshift(obj.data.maxwellCorrX,1),2),3);
    obj.data.maxwellCorrY =  ifftshift(ifftshift(ifftshift(obj.data.maxwellCorrY,1),2),3);
    obj.data.maxwellCorrZ =  ifftshift(ifftshift(ifftshift(obj.data.maxwellCorrZ,1),2),3);

    obj.data.sensMaps = ifftshift(ifftshift(ifftshift(obj.data.sensMaps,1),2),3);
    obj.data.x0 =  ifftshift(ifftshift(ifftshift(obj.data.x0,1),2),3);
end

function Out = reconGAMP(S)
        
        EstimIn = S.EstimIn;
        EstimOut = S.EstimOut;
        Op = S.Op;
        opt = S.opt;
        compute = S.compute;
        
        
        if strcmpi(compute, 'gpu')
            [opt, Op, EstimOut] = putOnGPU(opt,Op,EstimOut,[],[]);
        end
        
%         xhatGAMP = gampEst_memEff(EstimIn, EstimOut, Op, opt);
        xhatGAMP = gampEst(EstimIn, EstimOut, Op, opt);
        
        xhatGAMP = clearOpts(xhatGAMP);
        
        if strcmpi(compute, 'gpu')
            [opt, Op, EstimOut, ~, xhatGAMP] = putOnCPU(opt, Op, EstimOut, [], xhatGAMP);
        end
        
        Out.opt = opt;
        Out.Op = Op;
        Out.EstimOut = EstimOut;
        Out.xhatGAMP = xhatGAMP;

end

