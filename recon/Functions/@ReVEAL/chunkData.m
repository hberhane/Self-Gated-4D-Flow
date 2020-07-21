% =========================================================================
% CHUNKDATA Breaks the data in the image domain along the frequency encoded
% direction into several different overlapping pieces, or chunks, to speed
% up computation and allow large datasets to fit on GPU. 
% Call obj.chunkData(dim, numChunks); after importing the data and 
% estimation the sensitivity maps
% =========================================================================
function chunkData(obj, dim, nChunks, overlap_prc)

%% Define overlapping chunks
% *************************************************************************
% overlap_prc = 0.25; % fraction of overlap b/w adjacent chunks
dim = min(dim, length(size(obj.data.Yb)));

FE_size = size(obj.data.Yb, dim);
% chunkSize = (FE_size+nChunks-2) / (overlap_prc*(1-nChunks)+nChunks);
% chunkSize = (1-overlap_perc).*((FE_size./nChunks)+0.5*overlap_perc*L) + P1*((N./C)+0.5*P1*L);
% chunkSize = FE_size / (nChunks*(1-0.5*overlap_prc));
chunkSize = FE_size / (nChunks*(1-overlap_prc)+overlap_prc);


for i = 1 : nChunks
chunkInds(i,1) = (i-1)*chunkSize - (i-1)*overlap_prc*chunkSize + 1;
chunkInds(i,2) = min(i*chunkSize - (i-1)*(overlap_prc*chunkSize),FE_size);
end

% chunkInds(1,1) = 1;
% chunkInds(1,2) = chunkSize;
% for i = 2 : nChunks
%     chunkInds(i,1) = chunkInds(i-1,1) + (1-overlap_prc)*chunkSize - 1;
%     chunkInds(i,2) = min(chunkInds(i,1) + chunkSize, FE_size);
% end
chunkInds = round(chunkInds);
obj.data.chunkInds = chunkInds;
% *************************************************************************
    
if dim == 1
    % asymmetric echo
    asym_size = max(find(squeeze(sum(sum(sum(obj.data.sampB,4),3),2))==0));
    asym_percent = asym_size/FE_size;
end

for ch = 1 : nChunks
    
    if obj.options.is1Dir && obj.options.isPlanar
        
        % take data to image domain
        Yb = ifft2_shift(obj.data.Yb);
        Yx = ifft2_shift(obj.data.Yx);
        
        % Crop along frequency encoded direction
        % Data
        obj.data.chunk(ch).Yb = Yb(chunkInds(ch,1):chunkInds(ch,2),:,:,:,:);
        obj.data.chunk(ch).Yx = Yx(chunkInds(ch,1):chunkInds(ch,2),:,:,:,:);
        % Samples
        tmpSize = chunkInds(ch,2) - chunkInds(ch,1) + 1;
        obj.data.chunk(ch).sampB = obj.data.sampB(FE_size-tmpSize+1:FE_size,:,:,:);
        obj.data.chunk(ch).sampX = obj.data.sampX(FE_size-tmpSize+1:FE_size,:,:,:);
        % Weights
        if ~isempty(obj.data.weightsB)
            obj.data.chunk(ch).weightsB = obj.data.weightsB(FE_size-tmpSize+1:FE_size,:,:,:);
            obj.data.chunk(ch).weightsX = obj.data.weightsX(FE_size-tmpSize+1:FE_size,:,:,:);
        end
        % Maxwell correction maps
        obj.data.chunk(ch).maxwellCorrX = obj.data.maxwellCorrX(chunkInds(ch,1):chunkInds(ch,2),:,:);
        % Image initialization
        obj.data.chunk(ch).x0 = obj.data.x0(chunkInds(ch,1):chunkInds(ch,2),:,:,:,:);
        % Crop sensitivity maps if already estimated
        if ~isempty(obj.data.sensMaps)
            obj.data.chunk(ch).sensMaps = obj.data.sensMaps(chunkInds(ch,1):chunkInds(ch,2),:,:,:,:);
        end
        
        % convert back to k-space
        obj.data.chunk(ch).Yb = fft2_shift(obj.data.chunk(ch).Yb);
        obj.data.chunk(ch).Yx = fft2_shift(obj.data.chunk(ch).Yx);
        
        % reinforce asymmetric echo
        % samples
        obj.data.chunk(ch).sampB(1:round(tmpSize*asym_percent),:,:,:,:) = 0;
        obj.data.chunk(ch).sampX(1:round(tmpSize*asym_percent),:,:,:,:) = 0;
        % data
        obj.data.chunk(ch).Yb = bsxfun(@times,obj.data.chunk(ch).Yb,permute(obj.data.chunk(ch).sampB,[1,2,4,3]));
        obj.data.chunk(ch).Yx = bsxfun(@times,obj.data.chunk(ch).Yx,permute(obj.data.chunk(ch).sampX,[1,2,4,3]));
        % weights
        if ~isempty(obj.data.weightsB)
            obj.data.chunk(ch).weightsB = bsxfun(@times,obj.data.chunk(ch).weightsB,permute(obj.data.chunk(ch).sampB,[1,2,4,3]));
            obj.data.chunk(ch).weightsX = bsxfun(@times,obj.data.chunk(ch).weightsX,permute(obj.data.chunk(ch).sampX,[1,2,4,3]));
        end
        
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
        
        switch dim
            
            % Chunk along frequency encode dimension
            case 1
        
%                 % take data to image domain
                Yb = ifft3_shift(obj.data.Yb);
                Yx = ifft3_shift(obj.data.Yx);
                Yy = ifft3_shift(obj.data.Yy);
                Yz = ifft3_shift(obj.data.Yz);

%                 Yb = obj.data.Yb;
%                 Yx = obj.data.Yx;
%                 Yy = obj.data.Yy;
%                 Yz = obj.data.Yz;
                
                % Crop along frequency encoded direction
                % Data
                obj.data.chunk(ch).Yb = Yb(chunkInds(ch,1):chunkInds(ch,2),:,:,:,:);
                obj.data.chunk(ch).Yx = Yx(chunkInds(ch,1):chunkInds(ch,2),:,:,:,:);
                obj.data.chunk(ch).Yy = Yy(chunkInds(ch,1):chunkInds(ch,2),:,:,:,:);
                obj.data.chunk(ch).Yz = Yz(chunkInds(ch,1):chunkInds(ch,2),:,:,:,:);
                % Samples
                tmpSize = chunkInds(ch,2) - chunkInds(ch,1) + 1;
                obj.data.chunk(ch).sampB = obj.data.sampB(FE_size-tmpSize+1:FE_size,:,:,:);
                obj.data.chunk(ch).sampX = obj.data.sampX(FE_size-tmpSize+1:FE_size,:,:,:);
                obj.data.chunk(ch).sampY = obj.data.sampY(FE_size-tmpSize+1:FE_size,:,:,:);
                obj.data.chunk(ch).sampZ = obj.data.sampZ(FE_size-tmpSize+1:FE_size,:,:,:);
                % Weights
                if ~isempty(obj.data.weightsB)
                    obj.data.chunk(ch).weightsB = obj.data.weightsB(FE_size-tmpSize+1:FE_size,:,:,:);
                    obj.data.chunk(ch).weightsX = obj.data.weightsX(FE_size-tmpSize+1:FE_size,:,:,:);
                    obj.data.chunk(ch).weightsY = obj.data.weightsY(FE_size-tmpSize+1:FE_size,:,:,:);
                    obj.data.chunk(ch).weightsZ = obj.data.weightsZ(FE_size-tmpSize+1:FE_size,:,:,:);
                end
                % Maxwell correction maps
                obj.data.chunk(ch).maxwellCorrX = obj.data.maxwellCorrX(chunkInds(ch,1):chunkInds(ch,2),:,:);
                obj.data.chunk(ch).maxwellCorrY = obj.data.maxwellCorrY(chunkInds(ch,1):chunkInds(ch,2),:,:);
                obj.data.chunk(ch).maxwellCorrZ = obj.data.maxwellCorrZ(chunkInds(ch,1):chunkInds(ch,2),:,:);
                % Image initialization
                if ~isempty(obj.data.x0)
                    obj.data.chunk(ch).x0 = obj.data.x0(chunkInds(ch,1):chunkInds(ch,2),:,:,:,:);
                end
                % Crop sensitivity maps if already estimated
                if ~isempty(obj.data.sensMaps)
                    obj.data.chunk(ch).sensMaps = obj.data.sensMaps(chunkInds(ch,1):chunkInds(ch,2),:,:,:,:);
                end
%                 % convert back to k-space
                obj.data.chunk(ch).Yb = fft3_shift(obj.data.chunk(ch).Yb);
                obj.data.chunk(ch).Yx = fft3_shift(obj.data.chunk(ch).Yx);
                obj.data.chunk(ch).Yy = fft3_shift(obj.data.chunk(ch).Yy);
                obj.data.chunk(ch).Yz = fft3_shift(obj.data.chunk(ch).Yz);
                % reinforce asymmetric echo
                % samples
                obj.data.chunk(ch).sampB(1:round(tmpSize*asym_percent),:,:,:,:) = 0;
                obj.data.chunk(ch).sampX(1:round(tmpSize*asym_percent),:,:,:,:) = 0;
                obj.data.chunk(ch).sampY(1:round(tmpSize*asym_percent),:,:,:,:) = 0;
                obj.data.chunk(ch).sampZ(1:round(tmpSize*asym_percent),:,:,:,:) = 0;
                % data
                obj.data.chunk(ch).Yb = bsxfun(@times,obj.data.chunk(ch).Yb,permute(obj.data.chunk(ch).sampB,[1,2,3,5,4]));
                obj.data.chunk(ch).Yx = bsxfun(@times,obj.data.chunk(ch).Yx,permute(obj.data.chunk(ch).sampX,[1,2,3,5,4]));
                obj.data.chunk(ch).Yy = bsxfun(@times,obj.data.chunk(ch).Yy,permute(obj.data.chunk(ch).sampY,[1,2,3,5,4]));
                obj.data.chunk(ch).Yz = bsxfun(@times,obj.data.chunk(ch).Yz,permute(obj.data.chunk(ch).sampZ,[1,2,3,5,4]));
                
                % weights
                if ~isempty(obj.data.weightsB)
                    obj.data.chunk(ch).weightsB = bsxfun(@times,obj.data.chunk(ch).weightsB,obj.data.chunk(ch).sampB);
                    obj.data.chunk(ch).weightsX = bsxfun(@times,obj.data.chunk(ch).weightsX,obj.data.chunk(ch).sampX);
                    obj.data.chunk(ch).weightsY = bsxfun(@times,obj.data.chunk(ch).weightsY,obj.data.chunk(ch).sampY);
                    obj.data.chunk(ch).weightsZ = bsxfun(@times,obj.data.chunk(ch).weightsZ,obj.data.chunk(ch).sampZ);
                end
            
            % Chunk along phase encode dimension
            case 2
               
            % Chunk along partition encode dimension    
            case 3
            
            % Chunk along cardiac dimension    
            case 5
            
                % Crop along cardiac direction
                % Data
                obj.data.chunk(ch).Yb = obj.data.Yb(:,:,:,:,chunkInds(ch,1):chunkInds(ch,2));
                obj.data.chunk(ch).Yx = obj.data.Yx(:,:,:,:,chunkInds(ch,1):chunkInds(ch,2));
                obj.data.chunk(ch).Yy = obj.data.Yy(:,:,:,:,chunkInds(ch,1):chunkInds(ch,2));
                obj.data.chunk(ch).Yz = obj.data.Yz(:,:,:,:,chunkInds(ch,1):chunkInds(ch,2));
                % Samples
                obj.data.chunk(ch).sampB = obj.data.sampB(:,:,:,chunkInds(ch,1):chunkInds(ch,2));
                obj.data.chunk(ch).sampX = obj.data.sampX(:,:,:,chunkInds(ch,1):chunkInds(ch,2));
                obj.data.chunk(ch).sampY = obj.data.sampY(:,:,:,chunkInds(ch,1):chunkInds(ch,2));
                obj.data.chunk(ch).sampZ = obj.data.sampZ(:,:,:,chunkInds(ch,1):chunkInds(ch,2));
                % Weights
                if ~isempty(obj.data.weightsB)
                    obj.data.chunk(ch).weightsB = obj.data.weightsB(:,:,:,chunkInds(ch,1):chunkInds(ch,2));
                    obj.data.chunk(ch).weightsX = obj.data.weightsX(:,:,:,chunkInds(ch,1):chunkInds(ch,2));
                    obj.data.chunk(ch).weightsY = obj.data.weightsY(:,:,:,chunkInds(ch,1):chunkInds(ch,2));
                    obj.data.chunk(ch).weightsZ = obj.data.weightsZ(:,:,:,chunkInds(ch,1):chunkInds(ch,2));
                end
                % Maxwell correction maps
                obj.data.chunk(ch).maxwellCorrX = obj.data.maxwellCorrX(:,:,:);
                obj.data.chunk(ch).maxwellCorrY = obj.data.maxwellCorrY(:,:,:);
                obj.data.chunk(ch).maxwellCorrZ = obj.data.maxwellCorrZ(:,:,:);
                % Image initialization
                if ~isempty(obj.data.x0)
                    obj.data.chunk(ch).x0 = obj.data.x0(:,:,:,chunkInds(ch,1):chunkInds(ch,2));
                end
                % Crop sensitivity maps if already estimated
                if ~isempty(obj.data.sensMaps)
                    obj.data.chunk(ch).sensMaps = obj.data.sensMaps;
                end
                
        end
        
    end
end

end



