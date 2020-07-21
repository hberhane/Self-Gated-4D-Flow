function cropData(obj,crop_val)
%CROPDATA Crops the data in the image domain along the frequence encoded
%dimension to speed up computation.  Call obj.cropData(crop_val); after
%importing the data and estimation the sensitivty maps
%   Detailed explanation goes here

if obj.options.is1Dir && obj.options.isPlanar
    
    % take data to image domain
    obj.data.Yb = ifft2_shift(obj.data.Yb);
    obj.data.Yx = ifft2_shift(obj.data.Yx);

    % Crop along frequency encoded direction
    obj.data.Yb = obj.data.Yb(:,crop_val+1:end-crop_val,:,:);
    obj.data.Yx = obj.data.Yx(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampB = obj.data.sampB(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampX = obj.data.sampX(:,crop_val+1:end-crop_val,:,:);
    obj.data.maxwellCorrX = obj.data.maxwellCorrX(:,crop_val+1:end-crop_val,:,:);
    obj.data.x0 = obj.data.x0(:,crop_val+1:end-crop_val,:,:);
    
    % Crop sensitivity maps if already estimated
    if ~isempty(obj.data.sensMaps)
        obj.data.sensMaps = obj.data.sensMaps(:,crop_val+1:end-crop_val,:,:);
    end

    % convert back to k-space
    obj.data.Yb = fft2_shift(obj.data.Yb);
    obj.data.Yx = fft2_shift(obj.data.Yx);
    
elseif ~obj.options.is1Dir && obj.options.isPlanar
    % take data to image domain
    obj.data.Yb = ifft2_shift(obj.data.Yb);
    obj.data.Yx = ifft2_shift(obj.data.Yx);
    obj.data.Yy = ifft2_shift(obj.data.Yy);
    obj.data.Yz = ifft2_shift(obj.data.Yz);
    
    % Crop along frequency encoded direction
    obj.data.Yb = obj.data.Yb(:,crop_val+1:end-crop_val,:,:);
    obj.data.Yx = obj.data.Yx(:,crop_val+1:end-crop_val,:,:);
    obj.data.Yy = obj.data.Yy(:,crop_val+1:end-crop_val,:,:);
    obj.data.Yz = obj.data.Yz(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampB = obj.data.sampB(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampX = obj.data.sampX(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampY = obj.data.sampY(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampZ = obj.data.sampZ(:,crop_val+1:end-crop_val,:,:);
    obj.data.maxwellCorrX = obj.data.maxwellCorrX(:,crop_val+1:end-crop_val,:,:);
    obj.data.maxwellCorrY = obj.data.maxwellCorrY(:,crop_val+1:end-crop_val,:,:);
    obj.data.maxwellCorrZ = obj.data.maxwellCorrZ(:,crop_val+1:end-crop_val,:,:);
    obj.data.x0 = obj.data.x0(:,crop_val+1:end-crop_val,:,:);
    
    % Crop sensitivity maps if already estimated
    if ~isempty(obj.data.sensMaps)
        obj.data.sensMaps = obj.data.sensMaps(:,crop_val+1:end-crop_val,:,:);
    end

    % convert back to k-space
    obj.data.Yb = fft2_shift(obj.data.Yb);
    obj.data.Yx = fft2_shift(obj.data.Yx);
    obj.data.Yy = fft2_shift(obj.data.Yy);
    obj.data.Yz = fft2_shift(obj.data.Yz);
    
elseif ~obj.options.is1Dir && ~obj.options.isPlanar
    
    % Find asymmetric echo percent
    tmp = sum(sum(sum(obj.data.sampB, 2),3),4);
    asym_size = max(find(tmp==0));
    asym_percent = asym_size/size(obj.data.Yb,1);
    
    % take data to image domain
    obj.data.Yb = ifft3_shift(obj.data.Yb);
    obj.data.Yx = ifft3_shift(obj.data.Yx);
    obj.data.Yy = ifft3_shift(obj.data.Yy);
    obj.data.Yz = ifft3_shift(obj.data.Yz);
    
    crop_val_start = crop_val;
    crop_val_stop = crop_val;

    % Crop along frequency encoded direction
    obj.data.Yb = obj.data.Yb(crop_val_start+1:end-crop_val_stop,:,:,:,:);
    obj.data.Yx = obj.data.Yx(crop_val_start+1:end-crop_val_stop,:,:,:,:);
    obj.data.Yy = obj.data.Yy(crop_val_start+1:end-crop_val_stop,:,:,:,:);
    obj.data.Yz = obj.data.Yz(crop_val_start+1:end-crop_val_stop,:,:,:,:);
    obj.data.sampB = obj.data.sampB(crop_val_start+1:end-crop_val_stop,:,:,:);
    obj.data.sampX = obj.data.sampX(crop_val_start+1:end-crop_val_stop,:,:,:);
    obj.data.sampY = obj.data.sampY(crop_val_start+1:end-crop_val_stop,:,:,:);
    obj.data.sampZ = obj.data.sampZ(crop_val_start+1:end-crop_val_stop,:,:,:);
    obj.data.maxwellCorrX = obj.data.maxwellCorrX(crop_val_start+1:end-crop_val_stop,:,:,:);
    obj.data.maxwellCorrY = obj.data.maxwellCorrY(crop_val_start+1:end-crop_val_stop,:,:,:);
    obj.data.maxwellCorrZ = obj.data.maxwellCorrZ(crop_val_start+1:end-crop_val_stop,:,:,:);
    obj.data.x0 = obj.data.x0(crop_val_start+1:end-crop_val_stop,:,:,:);
    if ~isempty(obj.data.weightsB)
    obj.data.weightsB = obj.data.weightsB(crop_val_start+1:end-crop_val_stop,:,:,:);
    obj.data.weightsX = obj.data.weightsX(crop_val_start+1:end-crop_val_stop,:,:,:);
    obj.data.weightsY = obj.data.weightsY(crop_val_start+1:end-crop_val_stop,:,:,:);
    obj.data.weightsZ = obj.data.weightsZ(crop_val_start+1:end-crop_val_stop,:,:,:);
    end
    
    % Crop sensitivity maps if already estimated
    if ~isempty(obj.data.sensMaps)
        obj.data.sensMaps = obj.data.sensMaps(crop_val_start+1:end-crop_val_stop,:,:,:,:);
    end

    % convert back to k-space
    obj.data.Yb = fft3_shift(obj.data.Yb);
    obj.data.Yx = fft3_shift(obj.data.Yx);
    obj.data.Yy = fft3_shift(obj.data.Yy);
    obj.data.Yz = fft3_shift(obj.data.Yz);
    
    % reinforce asymmetric echo
    % Enforce Asymmetric Echo
    obj.data.sampB(1:round(size(obj.data.sampB,1)*asym_percent),:,:,:) = 0;
    obj.data.sampX(1:round(size(obj.data.sampX,1)*asym_percent),:,:,:) = 0;
    obj.data.sampY(1:round(size(obj.data.sampY,1)*asym_percent),:,:,:) = 0;
    obj.data.sampZ(1:round(size(obj.data.sampZ,1)*asym_percent),:,:,:) = 0;

    obj.data.Yb = bsxfun(@times,obj.data.Yb,permute(obj.data.sampB,[1,2,3,5,4]));
    obj.data.Yx = bsxfun(@times,obj.data.Yx,permute(obj.data.sampX,[1,2,3,5,4]));
    obj.data.Yy = bsxfun(@times,obj.data.Yy,permute(obj.data.sampY,[1,2,3,5,4]));
    obj.data.Yz = bsxfun(@times,obj.data.Yz,permute(obj.data.sampZ,[1,2,3,5,4]));

    obj.data.weightsB = bsxfun(@times,obj.data.weightsB,obj.data.sampB);
    obj.data.weightsX = bsxfun(@times,obj.data.weightsX,obj.data.sampX);
    obj.data.weightsY = bsxfun(@times,obj.data.weightsY,obj.data.sampY);
    obj.data.weightsZ = bsxfun(@times,obj.data.weightsZ,obj.data.sampZ);
end

end

