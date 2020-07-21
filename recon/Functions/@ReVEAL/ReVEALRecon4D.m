function ReVEALRecon4D(obj)

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
fprintf(sprintf('R = %s\n',num2str(R)))

% Check if multiple sets of sensitivity maps are used from espirit
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
if size(obj.data.x0, 5) == 1
obj.options.GAMPOptB.xhat0 = obj.data.x0(:);
obj.options.GAMPOptX.xhat0 = obj.data.x0(:);
obj.options.GAMPOptY.xhat0 = obj.data.x0(:);
obj.options.GAMPOptZ.xhat0 = obj.data.x0(:);
else
obj.options.GAMPOptB.xhat0 = obj.data.x0(:,:,:,:,1); obj.options.GAMPOptB.xhat0 = obj.options.GAMPOptB.xhat0(:);
obj.options.GAMPOptX.xhat0 = obj.data.x0(:,:,:,:,2); obj.options.GAMPOptX.xhat0 = obj.options.GAMPOptX.xhat0(:);
obj.options.GAMPOptY.xhat0 = obj.data.x0(:,:,:,:,3); obj.options.GAMPOptY.xhat0 = obj.options.GAMPOptY.xhat0(:);
obj.options.GAMPOptZ.xhat0 = obj.data.x0(:,:,:,:,4); obj.options.GAMPOptZ.xhat0 = obj.options.GAMPOptZ.xhat0(:);
end
    

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


% Reconstruct background image
fprintf('\npB...\n')
if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
    [obj.options.GAMPOptB,Op_b,EstimOut_b] = putOnGPU(obj.options.GAMPOptB,Op_b,EstimOut_b,[],[]);
end
if strcmp(obj.options.ReVEALOpts.data_size, 'large')
    xhatGAMP_b = gampEst_memEff(EstimIn,EstimOut_b,Op_b,obj.options.GAMPOptB);
else
    xhatGAMP_b = gampEst(EstimIn,EstimOut_b,Op_b,obj.options.GAMPOptB);
end
xhatGAMP_b = clearOpts(xhatGAMP_b);
if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
    [obj.options.GAMPOptB,Op_b,EstimOut_b,~,xhatGAMP_b] = putOnCPU(obj.options.GAMPOptB,Op_b,EstimOut_b,[],xhatGAMP_b);
end

% Reconstruct x encoded image
fprintf('\npX...\n')
if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
    [obj.options.GAMPOptX,Op_x,EstimOut_x] = putOnGPU(obj.options.GAMPOptX,Op_x,EstimOut_x,[],[]);
end
if strcmp(obj.options.ReVEALOpts.data_size, 'large')
    xhatGAMP_x = gampEst_memEff(EstimIn,EstimOut_x,Op_x,obj.options.GAMPOptX);
else
    xhatGAMP_x = gampEst(EstimIn,EstimOut_x,Op_x,obj.options.GAMPOptX);
end
xhatGAMP_x = clearOpts(xhatGAMP_x);
if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
    [obj.options.GAMPOptX,Op_x,EstimOut_x,~,xhatGAMP_x] = putOnCPU(obj.options.GAMPOptX,Op_x,EstimOut_x,[],xhatGAMP_x);
end

% Reconstruct y encoded image
fprintf('\npY...\n')
if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
    [obj.options.GAMPOptY,Op_y,EstimOut_y] = putOnGPU(obj.options.GAMPOptY,Op_y,EstimOut_y,[],[]);
end
if strcmp(obj.options.ReVEALOpts.data_size, 'large')
    xhatGAMP_y = gampEst_memEff(EstimIn,EstimOut_y,Op_y,obj.options.GAMPOptY);
else
    xhatGAMP_y = gampEst(EstimIn,EstimOut_y,Op_y,obj.options.GAMPOptY);
end
xhatGAMP_y = clearOpts(xhatGAMP_y);
if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
    [obj.options.GAMPOptY,Op_y,EstimOut_y,~,xhatGAMP_y] = putOnCPU(obj.options.GAMPOptY,Op_y,EstimOut_y,[],xhatGAMP_y);
end

% Reconstruct z encoded image
fprintf('\npZ...\n')
if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
    [obj.options.GAMPOptZ,Op_z,EstimOut_z] = putOnGPU(obj.options.GAMPOptZ,Op_z,EstimOut_z,[],[]);
end
if strcmp(obj.options.ReVEALOpts.data_size, 'large')
    xhatGAMP_z = gampEst_memEff(EstimIn,EstimOut_z,Op_z,obj.options.GAMPOptZ);
else
    xhatGAMP_z = gampEst(EstimIn,EstimOut_z,Op_z,obj.options.GAMPOptZ);
end
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


% apply fftshifts and restimate sensitivity maps after BPO correction
applyIfftshift(obj);
obj.estimateSensMaps(1); % averaging all encodings
obj.data.sensMaps = obj.data.sensMaps*sqrt(R);
applyFftshift(obj);

% Re-estimate sigmaSq from initial recon
% avg_image_b = reshape(xhatGAMP_b.xhat,size(obj.data.sampB));
% avg_image_x = reshape(xhatGAMP_x.xhat,size(obj.data.sampX));
% avg_image_y = reshape(xhatGAMP_y.xhat,size(obj.data.sampY));
% avg_image_z = reshape(xhatGAMP_z.xhat,size(obj.data.sampZ));
% 
% % avg_image = cat(4,avg_image_b,avg_image_x,avg_image_y,avg_image_z);
% % avg_image = mean(avg_image, 4);
% % avg_image = avg_image(:);
% % avg_image(avg_image<0.01*max(avg_image)) = [];
% % sigma = median(abs(avg_image));
% % tau = 0.15; % normally 0.15
% % sigmaSq = (sigma*tau)^2;
% % obj.options.ReVEALOpts.sigmaSq = sigmaSq;
% % obj.options.ReVEALOpts.sigmaSq = 4*obj.options.ReVEALOpts.wvar;
% 
% avg_mag = squeeze(mean(abs(cat(5,avg_image_b,avg_image_x,avg_image_y,avg_image_z)), 4));
% avg_mag(avg_mag < 0.02*max(avg_mag(:))) = 0;
% dif(:,:,:,1) = avg_mag(:,:,:,2) - avg_mag(:,:,:,1);
% dif(:,:,:,2) = avg_mag(:,:,:,3) - avg_mag(:,:,:,1);
% dif(:,:,:,3) = avg_mag(:,:,:,4) - avg_mag(:,:,:,1);
% dif = cat(1, dif(:,:,:,1), dif(:,:,:,2), dif(:,:,:,3));
% dif(dif == 0) = [];

% ********************************
% OPTION 1
% sigsq = 2*(std(dif(:))^2);
% ********************************

% ********************************
% OPTION 2 
% factor = 3; % user-defined, 3 seems to work
% sigsq = (factor*prctile(dif(:), 95))^2; % aim for around sigmasq ~ 0.01
% ********************************
% sigsq = prctile(dif, 90);
% sig(sig == 0) = [];

% hardcode sigma_sq
% obj.options.ReVEALOpts.sigmaSq = 0.01;



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
    
    % ---------------------------------------------------------------------
    % After so many iterations, identify and downweight readout
    % lines which do not conform well to the fit
    if ind == 4
        % Compute residual of the data fidelty term
        xhat_b = reshape(xhatGAMP_b.xhat, [size(obj.data.sampB)]);
        xhat_x = reshape(xhatGAMP_x.xhat, [size(obj.data.sampX)]);
        xhat_y = reshape(xhatGAMP_y.xhat, [size(obj.data.sampY)]);
        xhat_z = reshape(xhatGAMP_z.xhat, [size(obj.data.sampZ)]);
        
        Ax_b = bsxfun(@times, permute(obj.data.weightsB, [1,2,3,5,4]), bsxfun(@times, permute(obj.data.sampB, [1,2,3,5,4]), fft3_shift(bsxfun(@times, obj.data.sensMaps, permute(xhat_b, [1,2,3,5,4])))));
        Ax_x = bsxfun(@times, permute(obj.data.weightsX, [1,2,3,5,4]), bsxfun(@times, permute(obj.data.sampX, [1,2,3,5,4]), fft3_shift(bsxfun(@times, obj.data.sensMaps, permute(xhat_x, [1,2,3,5,4])))));
        Ax_y = bsxfun(@times, permute(obj.data.weightsY, [1,2,3,5,4]), bsxfun(@times, permute(obj.data.sampY, [1,2,3,5,4]), fft3_shift(bsxfun(@times, obj.data.sensMaps, permute(xhat_y, [1,2,3,5,4])))));
        Ax_z = bsxfun(@times, permute(obj.data.weightsZ, [1,2,3,5,4]), bsxfun(@times, permute(obj.data.sampZ, [1,2,3,5,4]), fft3_shift(bsxfun(@times, obj.data.sensMaps, permute(xhat_z, [1,2,3,5,4])))));
        
        wResidual = cat(6, Ax_b-obj.data.Yb, Ax_x-obj.data.Yx, Ax_y-obj.data.Yy, Ax_z-obj.data.Yz);
        tmp = abs(wResidual);
        tmp(tmp == 0) = [];
                    
        % zscore to identify any outliers
        zscore = (abs(wResidual) - mean(tmp))/std(tmp);
        thresh = 5;
        keepIndx = zeros(size(wResidual),'logical');
        keepIndx(abs(zscore) < thresh) = 1;
        keepIndx = squeeze(prod(keepIndx, 4));
        
        % reinforce sampling pattern
        keepIndx_b = bsxfun(@times, obj.data.sampB, keepIndx(:,:,:,:,1));
        keepIndx_x = bsxfun(@times, obj.data.sampX, keepIndx(:,:,:,:,2));
        keepIndx_y = bsxfun(@times, obj.data.sampY, keepIndx(:,:,:,:,3));
        keepIndx_z = bsxfun(@times, obj.data.sampZ, keepIndx(:,:,:,:,4));
        
        sampsKept = sum(keepIndx_b(:))+sum(keepIndx_x(:))+sum(keepIndx_y(:))+sum(keepIndx_z(:));
        sampsOrig = sum(obj.data.sampB(:))+sum(obj.data.sampX(:))+sum(obj.data.sampY(:))+sum(obj.data.sampZ(:));
        removePerc = (1 - sampsKept/sampsOrig)*100;
        disp(['NOT IMPLEMENTED: FOR TESTING PURPOSES CURRENTLY']);
        disp([num2str(removePerc),'% of samples downweighted at ',num2str(ind),'th iteration']);
        
        % ****************************************************************
        
        
%         % set weights of non-conforming samples to zero
%         obj.data.weightsB = bsxfun(@times, keepIndx_b, obj.data.weightsB);
%         obj.data.weightsX = bsxfun(@times, keepIndx_x, obj.data.weightsX);
%         obj.data.weightsY = bsxfun(@times, keepIndx_y, obj.data.weightsY);
%         obj.data.weightsZ = bsxfun(@times, keepIndx_z, obj.data.weightsZ);
%         
%         % reweight data
%         obj.data.Yb = bsxfun(@times,permute(keepIndx_b,[1,2,3,5,4]) ,obj.data.Yb);
%         obj.data.Yx = bsxfun(@times,permute(keepIndx_b,[1,2,3,5,4]) ,obj.data.Yx);
%         obj.data.Yy = bsxfun(@times,permute(keepIndx_b,[1,2,3,5,4]) ,obj.data.Yy);
%         obj.data.Yz = bsxfun(@times,permute(keepIndx_b,[1,2,3,5,4]) ,obj.data.Yz);
%         
%         % Downsample Data
%         dataB = downsample_data(obj.data.Yb,obj.data.sampB);
%         dataX = downsample_data(obj.data.Yx,obj.data.sampX);
%         dataY = downsample_data(obj.data.Yy,obj.data.sampY);
%         dataZ = downsample_data(obj.data.Yz,obj.data.sampZ);
%         
%         % Redefine operators
%         % Create new Input Estimation Class with Maxwell Corrections
%         pMRI_b = pMRI_Op_3D_t(obj.data.sensMaps,obj.data.sampB,'uniform_var',...
%             obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute...
%             ,'precision',obj.options.ReVEALOpts.precision,'weights',obj.data.weightsB);
%         pMRI_x = pMRI_Op_3D_t(bsxfun(@times,obj.data.sensMaps,obj.data.maxwellCorrX),obj.data.sampX,'uniform_var',...
%             obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute...
%             ,'precision',obj.options.ReVEALOpts.precision,'weights',obj.data.weightsX);
%         pMRI_y = pMRI_Op_3D_t(bsxfun(@times,obj.data.sensMaps,obj.data.maxwellCorrY),obj.data.sampY,'uniform_var',...
%             obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute...
%             ,'precision',obj.options.ReVEALOpts.precision,'weights',obj.data.weightsY);
%         pMRI_z = pMRI_Op_3D_t(bsxfun(@times,obj.data.sensMaps,obj.data.maxwellCorrZ),obj.data.sampZ,'uniform_var',...
%             obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute...
%             ,'precision',obj.options.ReVEALOpts.precision,'weights',obj.data.weightsZ);
%         
%         % Create nd-DWT Linear Transform
%         W = ndDWTLinTrans4D(obj.options.SparseTrans,size(obj.data.sampB),'uniform_var',...
%             obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute,...
%             'precision',obj.options.ReVEALOpts.precision);
% 
%         % Concatonate All The Linear Transform Operator Together
%         Op_b = LinTransConcat({pMRI_b;W},[1,1],obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
%         Op_x = LinTransConcat({pMRI_x;W},[1,1],obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
%         Op_y = LinTransConcat({pMRI_y;W},[1,1],obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
%         Op_z = LinTransConcat({pMRI_z;W},[1,1],obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
%         
%         % Create Output Estimation classes
%         MeasEstimOut_b = CAwgnEstimOut(dataB,obj.options.ReVEALOpts.wvar,obj.options.ReVEALOpts.MAP);
%         MeasEstimOut_x = CAwgnEstimOut(dataX,obj.options.ReVEALOpts.wvar,obj.options.ReVEALOpts.MAP);
%         MeasEstimOut_y = CAwgnEstimOut(dataY,obj.options.ReVEALOpts.wvar,obj.options.ReVEALOpts.MAP);
%         MeasEstimOut_z = CAwgnEstimOut(dataZ,obj.options.ReVEALOpts.wvar,obj.options.ReVEALOpts.MAP);
%         
%         % Set the Regularization per band
%         lambda = setLambda(size(obj.data.sampB),obj.options.ReVEALOpts.lambda0);
%         AnaEstimOut1 = CplxLaplaceEstimOut(lambda);
%         
%         % Create Output Estimation class
%         EstimOut_b = EstimOutConcat({MeasEstimOut_b;AnaEstimOut1},[pMRI_b.M,W.M],...
%             obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
%         EstimOut_x = EstimOutConcat({MeasEstimOut_x;AnaEstimOut1},[pMRI_x.M,W.M],...
%             obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
%         EstimOut_y = EstimOutConcat({MeasEstimOut_y;AnaEstimOut1},[pMRI_y.M,W.M],...
%             obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
%         EstimOut_z = EstimOutConcat({MeasEstimOut_z;AnaEstimOut1},[pMRI_z.M,W.M],...
%             obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
%         
%         % Clear some memory up
%         clear meanB meanX meanY meanZ varB varX varY varZ pvarb pvarx pvary pvarz...
%             MeasEstimOut_b MeasEstimOut_x MeasEstimOut_y MeasEstimOut_z...
%             AnaEstimOut1 pMRI_b pMRI_x pMRI_y pMRI_z W dataB dataX dataY dataZ ...
%             lambda;
% %         % ********************************************************************
    end
    % ----------------------------------------------------------------------------------------------
       
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
    if strcmp(obj.options.ReVEALOpts.data_size, 'large')
        xhatGAMP_b = gampEst_memEff(EstimIn_b,EstimOut_b,Op_b,obj.options.GAMPOptB);
    else
        xhatGAMP_b = gampEst(EstimIn_b,EstimOut_b,Op_b,obj.options.GAMPOptB);
    end
   
    xhatGAMP_b = clearOpts(xhatGAMP_b);
    if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
        [obj.options.GAMPOptB,Op_b,EstimOut_b,~,xhatGAMP_b] = putOnCPU(obj.options.GAMPOptB,Op_b,EstimOut_b,EstimIn_b,xhatGAMP_b);
    end

    % Reconstruct x encoded image
    fprintf('\npX...\n')
    if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
        [obj.options.GAMPOptX,Op_x,EstimOut_x] = putOnGPU(obj.options.GAMPOptX,Op_x,EstimOut_x,EstimIn_x,[]);
    end
    if strcmp(obj.options.ReVEALOpts.data_size, 'large')
        xhatGAMP_x = gampEst_memEff(EstimIn_x,EstimOut_x,Op_x,obj.options.GAMPOptX);
    else
        xhatGAMP_x = gampEst(EstimIn_x,EstimOut_x,Op_x,obj.options.GAMPOptX);
    end
    xhatGAMP_x = clearOpts(xhatGAMP_x);
    if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
        [obj.options.GAMPOptX,Op_x,EstimOut_x,~,xhatGAMP_x] = putOnCPU(obj.options.GAMPOptX,Op_x,EstimOut_x,EstimIn_x,xhatGAMP_x);
    end

    % Reconstruct y encoded image
    fprintf('\npY...\n')
    if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
        [obj.options.GAMPOptY,Op_y,EstimOut_y] = putOnGPU(obj.options.GAMPOptY,Op_y,EstimOut_y,EstimIn_y,[]);
    end
    if strcmp(obj.options.ReVEALOpts.data_size, 'large')
        xhatGAMP_y = gampEst_memEff(EstimIn_y,EstimOut_y,Op_y,obj.options.GAMPOptY);
    else
        xhatGAMP_y = gampEst(EstimIn_y,EstimOut_y,Op_y,obj.options.GAMPOptY);
    end
    xhatGAMP_y = clearOpts(xhatGAMP_y);
    if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
        [obj.options.GAMPOptY,Op_y,EstimOut_y,~,xhatGAMP_y] = putOnCPU(obj.options.GAMPOptY,Op_y,EstimOut_y,EstimIn_y,xhatGAMP_y);
    end

    % Reconstruct z encoded image
    fprintf('\npZ...\n')
    if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
        [obj.options.GAMPOptZ,Op_z,EstimOut_z] = putOnGPU(obj.options.GAMPOptZ,Op_z,EstimOut_z,EstimIn_z,[]);
    end
    if strcmp(obj.options.ReVEALOpts.data_size, 'large')
        xhatGAMP_z = gampEst_memEff(EstimIn_z,EstimOut_z,Op_z,obj.options.GAMPOptZ);
    else
        xhatGAMP_z = gampEst(EstimIn_z,EstimOut_z,Op_z,obj.options.GAMPOptZ);
    end
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

obj.outputs.sensMaps = obj.data.sensMaps/sqrt(R);

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
%    warning('Change to ifftshift')
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
    
    
    saveDir_full = [obj.saveDir,'\Iter\',name];
    if ~exist(saveDir_full,'dir')
        mkdir(saveDir_full);
    end
        
    obj.saveGIFs(fullfile(saveDir_full, name))
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