function estimateSensMaps3D(obj,avg_all)
%ESTIMATESENSITIVITYMAPS Summary of this function goes here
%   Detailed explanation goes here

% fprintf('\n')
% star_display('Estimating Sensitvity Maps',0)
% tic;

if nargin < 2
    avg_all = 0;
end

% *** UPDATED averaging of all maps ********
% (1) Time-averaged images for each encoding
% Background
if ~isempty(obj.data.weightsB)
    weights = obj.data.weightsB;
else
    weights = obj.data.sampB;
end
sum_weights = sum(weights.^1,4);
sum_weights(sum_weights==0) = inf;
weights_sca = bsxfun(@rdivide,weights.^1,sum_weights);
avg_imageb = ifft3_shift(sum(bsxfun(@times,obj.data.Yb,permute(weights_sca,[1,2,3,5,4])),5));
% X encoding
if ~isempty(obj.data.weightsB)
    weights = obj.data.weightsX;
else
    weights = obj.data.sampX;
end
sum_weights = sum(weights.^1,4);
sum_weights(sum_weights==0) = inf;
weights_sca = bsxfun(@rdivide,weights.^1,sum_weights);
avg_imagex = ifft3_shift(sum(bsxfun(@times,obj.data.Yx,permute(weights_sca,[1,2,3,5,4])),5));
% Y encoding
if ~isempty(obj.data.weightsB)
    weights = obj.data.weightsY;
else
    weights = obj.data.sampY;
end
sum_weights = sum(weights.^1,4);
sum_weights(sum_weights==0) = inf;
weights_sca = bsxfun(@rdivide,weights.^1,sum_weights);
avg_imagey = ifft3_shift(sum(bsxfun(@times,obj.data.Yy,permute(weights_sca,[1,2,3,5,4])),5));
% Z encoding
if ~isempty(obj.data.weightsB)
    weights = obj.data.weightsZ;
else
    weights = obj.data.sampZ;
end
sum_weights = sum(weights.^1,4);
sum_weights(sum_weights==0) = inf;
weights_sca = bsxfun(@rdivide,weights.^1,sum_weights);
avg_imagez = ifft3_shift(sum(bsxfun(@times,obj.data.Yz,permute(weights_sca,[1,2,3,5,4])),5));
% (2) Average images for initial guess
avg_image = (avg_imageb + avg_imagex + avg_imagey + avg_imagez) / 4;
% (3) Sensitivity map estimation
p.mthd = 3; % '1' for espirit, '2' time-average spirit, '3' for walsh
p.reEst = 0; % Res-estimating sensitivities
p.fil = 3;
p.opt = 2;
% Use all 4 of the maps or just the background
if avg_all
    [senseB,~] = WalshCoilCombine3D(avg_imageb,p);
    [senseX,~] = WalshCoilCombine3D(avg_imagex,p);
    [senseY,~] = WalshCoilCombine3D(avg_imagey,p);
    [senseZ,~] = WalshCoilCombine3D(avg_imagez,p);
    % Remove background phase from maps
    senseX = bsxfun(@times, senseX, conj(obj.data.maxwellCorrX));
    senseY = bsxfun(@times, senseY, conj(obj.data.maxwellCorrY));
    senseZ = bsxfun(@times, senseZ, conj(obj.data.maxwellCorrZ));
    % Average maps
    obj.data.sensMaps = (senseB + senseX + senseY + senseZ) / 4;
else
    [obj.data.sensMaps,~] = WalshCoilCombine3D(avg_image, p);
end



% grow = avg_pattern(:,round(size(avg_pattern,2)/2),round(size(avg_pattern,3)/2));
% grow(grow==Inf)=NaN;
% dim_range1 = grow_dim(grow);

% grow = avg_pattern(round(size(avg_pattern,1)/2),:,round(size(avg_pattern,3)/2));
% grow(grow==Inf)=NaN;
% dim_range2 = grow_dim(grow.')

% grow = squeeze(avg_pattern(dim_range1(1),round(size(avg_pattern,2)/2),:));
% grow(grow==Inf)=NaN;
% dim_range3 = grow_dim(grow);

% avg_kspace = fft3_shift(avg_image);
% window = hamming(dim_range1(2)-dim_range1(1)+1);
% avg_kspace(dim_range1(1):dim_range1(2),:,:,:) = bsxfun(@times,avg_kspace(dim_range1(1):dim_range1(2),:,:,:),window);
% avg_image = ifft3_shift(avg_kspace);

% Initialize with SOS time averaged image
% avg_imageb = ifft3_shift(obj.data.Yb);
% avg_imagex = ifft3_shift(obj.data.Yx);
% avg_imagey = ifft3_shift(obj.data.Yy);
% avg_imagez = ifft3_shift(obj.data.Yz);

% avg_image = sum(avg_imageb,5)+sum(avg_imagex,5)+sum(avg_imagey,5)+sum(avg_imagez,5);
% avg_image = avg_image/(4*size(obj.data.Yb,5));

% SoS of time-averaged (4 encodings) images for initial guess
obj.data.x0 = sqrt(sum(abs(avg_image).^2,4));
obj.data.x0 = repmat(obj.data.x0,[1,1,1,size(obj.data.sampB,4)]);
% 
% % Estimate sensitivity maps
% p.mthd = 3; % '1' for espirit, '2' time-average spirit, '3' for walsh
% p.reEst = 0; % Res-estimating sensitivities
% p.fil = 3;
% p.opt = 2;
% [obj.data.sensMaps,~] = WalshCoilCombine3D(avg_image,p);

% *******************************************
% [senseB,~] = WalshCoilCombine3D(avg_imageb,p);
% [senseX,~] = WalshCoilCombine3D(avg_imagex,p);
% [senseY,~] = WalshCoilCombine3D(avg_imagey,p);
% [senseZ,~] = WalshCoilCombine3D(avg_imagez,p);
% 
% senseAvg_maps = (senseB + bsxfun(@times, senseX, conj(obj.data.maxwellCorrX)) ...
%                         + bsxfun(@times, senseY, conj(obj.data.maxwellCorrY)) ...
%                         + bsxfun(@times, senseZ, conj(obj.data.maxwellCorrZ))) / 4;
% coil = 1;
% figure
% subplot(231); imagesc(squeeze(abs(senseB(:,:,39,coil))),[0,1]); axis('off','image'); title('b');
% subplot(232); imagesc(squeeze(abs(senseX(:,:,39,coil))),[0,1]); axis('off','image'); title('x');
% subplot(233); imagesc(squeeze(abs(senseY(:,:,39,coil))),[0,1]); axis('off','image'); title('y');
% subplot(234); imagesc(squeeze(abs(senseZ(:,:,39,coil))),[0,1]); axis('off','image'); title('z');
% subplot(235); imagesc(squeeze(abs(senseAvg_kspace(:,:,39,coil))),[0,1]); axis('off','image'); title('Avg kspace');
% subplot(236); imagesc(squeeze(abs(senseAvg_maps(:,:,39,coil))),[0,1]); axis('off','image'); title('Avg maps');
% suptitle('magnitude');
% 
% 
% figure
% subplot(231); imagesc(squeeze(angle(senseB(:,:,39,coil)))); axis('off','image'); title('b');
% subplot(232); imagesc(squeeze(angle(senseX(:,:,39,coil)))); axis('off','image'); title('x');
% subplot(233); imagesc(squeeze(angle(senseY(:,:,39,coil)))); axis('off','image'); title('y');
% subplot(234); imagesc(squeeze(angle(senseZ(:,:,39,coil)))); axis('off','image'); title('z');
% subplot(235); imagesc(squeeze(angle(senseAvg_kspace(:,:,39,coil)))); axis('off','image'); title('Avg kspace');
% subplot(236); imagesc(squeeze(angle(senseAvg_maps(:,:,39,coil)))); axis('off','image'); title('Avg maps');
% suptitle('phase');
% 
% figure
% subplot(231); imagesc(squeeze(abs(senseB(:,:,39,coil))),[0,1]); axis('off','image'); title('b');
% subplot(232); imagesc(squeeze(abs(senseX(:,:,39,coil)-senseB(:,:,39,coil))),[0,1]); axis('off','image'); title('x');
% subplot(233); imagesc(squeeze(abs(senseY(:,:,39,coil)-senseB(:,:,39,coil))),[0,1]); axis('off','image'); title('y');
% subplot(234); imagesc(squeeze(abs(senseZ(:,:,39,coil)-senseB(:,:,39,coil))),[0,1]); axis('off','image'); title('z');
% subplot(235); imagesc(squeeze(abs(senseAvg_kspace(:,:,39,coil)-senseB(:,:,39,coil))),[0,1]); axis('off','image'); title('Avg kspace');
% subplot(236); imagesc(squeeze(abs(senseAvg_maps(:,:,39,coil)-senseB(:,:,39,coil))),[0,1]); axis('off','image'); title('Avg maps');
% suptitle('dif mag wrt B')
%                     
% figure
% subplot(231); imagesc(squeeze(angle(senseB(:,:,39,coil)))); axis('off','image'); title('b');
% subplot(232); imagesc(squeeze(angle(senseX(:,:,39,coil)-senseB(:,:,39,coil)))); axis('off','image'); title('x');
% subplot(233); imagesc(squeeze(angle(senseY(:,:,39,coil)-senseB(:,:,39,coil)))); axis('off','image'); title('y');
% subplot(234); imagesc(squeeze(angle(senseZ(:,:,39,coil)-senseB(:,:,39,coil)))); axis('off','image'); title('z');
% subplot(235); imagesc(squeeze(angle(senseAvg_kspace(:,:,39,coil)-senseB(:,:,39,coil)))); axis('off','image'); title('Avg kspace');
% subplot(236); imagesc(squeeze(angle(senseAvg_maps(:,:,39,coil)-senseB(:,:,39,coil)))); axis('off','image'); title('Avg maps');
% suptitle('dif phase wrt B')                    
                    
                    
% *********************************************                    
                    

% [obj.data.sensMaps,~] = coilSen3D(obj.data.Yb,obj.data.sampB,p);
% sensX = coilSen3D(obj.data.Yx,obj.data.sampX,p);
% sensY = coilSen3D(obj.data.Yy,obj.data.sampY,p);
% sensZ = coilSen3D(obj.data.Yz,obj.data.sampZ,p);
% obj.data.sensMaps = (obj.data.sensMaps + sensX + sensY +sensZ)/4;

% p.mthd   = 2; %'1' espirit, '2' time-average espirit, '3' walsh
% p.fil    = [6,6,6];  % size of kernal for eSPIRiT or size of filter for 'Walsh'; use [6,6,6] for eSPIRiT and [3,1,1] for Walsh
% p.ACSsz  = [24,24,24]; % size of the k-space block used for eSPIRiT
% p.eSRT   = 0.9;
% p.eSmaps = 1; %number of Espirit sensitivity maps
% p.avgPhs = 1; % Assign the phase to time-average image. 1: yes, 0: no
% p.ACSco  = [1/sqrt(2),1/sqrt(2)]; % Size of the cut-off filter used in Walsh ???????????????? use [1/2, 1/2 ]
% p.espirt = 0.5; %
% [obj.data.sensMaps,~]  = coilSen3D(obj.data.Yb, obj.data.sampB, p);
% 
% if size(obj.data.sensMaps, 5) == 2
%     
%     S1_tmp = obj.data.sensMaps(:,:,:,:,1); S1_tmpLin = S1_tmp(:);
%     S2_tmp = obj.data.sensMaps(:,:,:,:,2); S2_tmpLin = S2_tmp(:);
%     
%     ratio = norm(S1_tmpLin, 2) / norm(S2_tmpLin, 2);
%     
%     if ratio < 1
%         
%         obj.data.sensMaps(:,:,:,:,1) = S2_tmp;
%         obj.data.sensMaps(:,:,:,:,2) = S1_tmp;
%         
%     end
% end
        


% p.mthd   = 2; %'1' espirit, '2' Walsh
% p.fil    = [3,3,1];  % size of kernal for eSPIRiT or size of filter for 'Walsh'; use [6,6,6] for eSPIRiT and [3,1,1] for Walsh
% p.ACSsz  = [128,80,6]; % size of the k-space block used for eSPIRiT
% p.eSRT   = 0.95;
% p.eSmaps = 1; %number of Espirit sensitivity maps
% p.avgPhs = 1; % Assign the phase to time-average image. 1: yes, 0: no
% p.ACSco  = [1/sqrt(2),1/sqrt(2)]; % Size of the cut-off filter used in Walsh ???????????????? use [1/2, 1/2 ]
% [obj.data.sensMaps,obj.data.x0]  = coilSen(obj.data.Yb, obj.data.sampB, p);
% obj.data.sensMaps = obj.data.sensMaps(:,:,:,:,2);











% % I guess write out the code here b/c obj.estimateSensMaps() is not written
% % in such a way to do this properly (integrate later)
% %(1) Images after 1st iteration
% im_b = reshape(xhatGAMP_b.xhat,size(obj.data.sampB));
% im_x = reshape(xhatGAMP_x.xhat,size(obj.data.sampX));
% im_y = reshape(xhatGAMP_y.xhat,size(obj.data.sampY));
% im_z = reshape(xhatGAMP_z.xhat,size(obj.data.sampZ));
% % (2) Apply background phase correction
% im_b = im_b;
% im_x = bsxfun(@times, im_x, conj(obj.data.maxwellCorrX));
% im_y = bsxfun(@times, im_y, conj(obj.data.maxwellCorrY));
% im_z = bsxfun(@times, im_z, conj(obj.data.maxwellCorrZ));
% % (3) Multiply by current sensitivity maps
% im_b = bsxfun(@times, im_b, permute(obj.data.sensMaps, [1,2,3,5,4]));
% im_x = bsxfun(@times, im_x, permute(obj.data.sensMaps, [1,2,3,5,4]));
% im_y = bsxfun(@times, im_y, permute(obj.data.sensMaps, [1,2,3,5,4]));
% im_z = bsxfun(@times, im_z, permute(obj.data.sensMaps, [1,2,3,5,4]));
% % (4) Take to k-space
% k_b = fft3_shift(im_b);
% k_x = fft3_shift(im_x);
% k_y = fft3_shift(im_y);
% k_z = fft3_shift(im_z);
% % (5) Concatentate and average the kspace
% k = cat(4, k_b, k_x, k_y, k_z);
% k = squeeze(mean(k, 4));
% % (6) Crop kspace (trust center points more)
% N = size(k); M = round(size(k)/4);
% k = k((N(1)/2+1)-M(1):(N(1)/2+1)+M(1), (N(2)/2+1)-M(2):(N(2)/2+1)+M(2), (N(3)/2+1)-M(3):(N(3)/2+1)+M(3), :); 
% 
% % % k(:, [1:floor(size(k,2)/4),ceil(3*(size(k,2)/4)):end], [1:floor(size(k,3)/4),ceil(3*(size(k,3)/4)):end]) = 0;
% % k(:, [1:floor(size(k,2)/4),ceil(3*(size(k,2)/4)):end], :, :) = 0;
% % k(:, :, [1:floor(size(k,3)/4),ceil(3*(size(k,3)/4)):end], :) = 0;
% % (7) Apply windowing in kspace
% % ham2(1,:,1) = hamming(ceil(3*(size(k,2)/4))-floor(size(k,2)/4));
% % ham3(1,1,:) = hamming(ceil(3*(size(k,3)/4))-floor(size(k,3)/4));
% ham1(:,1,1) = hamming(size(k, 1));
% ham2(1,:,1) = hamming(size(k, 2));
% ham3(1,1,:) = hamming(size(k, 3));
% window = bsxfun(@times, ham1, bsxfun(@times, ham2, ham3));
% 
% % window = padarray(window, [0, floor(size(k,2)/4), floor(size(k,3)/4)], 'both');
% k = bsxfun(@times, k, window);
% M = size(k);
% k = padarray(k, ceil((N-M)/2),  'pre');
% k = padarray(k, floor((N-M)/2), 'post');
% % (7) Bring back to image domain
% im = ifft3_shift(k);
% % (7) restimate sensitivity maps
% p.mthd = 3; % '1' for espirit, '2' time-average spirit, '3' for walsh
% p.reEst = 0; % Res-estimating sensitivities
% p.fil = 3;
% p.opt = 2;
% [obj.data.sensMaps,~] = WalshCoilCombine3D(im,p);
% obj.data.sensMaps = obj.data.sensMaps*sqrt(R);




end

function dim_range = grow_dim(grow_1)
    center = round(size(grow_1,1)/2)+1;
    col_ind = 0;
    for ind = 1:size(grow_1,1)/2
        col = grow_1(center-ind:center+ind-1);
        if ~isnan(sum(col))
            col_ind = ind;
        end
    end

    dim_range = [center-col_ind, min(length(grow_1),center+col_ind)];
end
