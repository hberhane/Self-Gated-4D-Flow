function [ResCol,ResLin,ResPar,Rate_s,Rate_w,NImageCols,NImageLins,NImagePars] = constructArray_v2(FlowSG4D_Outputs, opt, name, saveLocation, param, raw_header)


SG_Reveal4D_Data = FlowSG4D_Outputs;
clear FlowSG4D_Outputs


ar_size = SG_Reveal4D_Data.ArraySize;
ar_ind = SG_Reveal4D_Data.Indices;
ar_ind = [ar_ind,(1:size(ar_ind,1)).'];
data = SG_Reveal4D_Data.Data;
rWeights = SG_Reveal4D_Data.rWeights;
rWeights = squeeze(single(rWeights));
    
%     ar_size(end) = 1;

crop = max([round(ar_size(1) / 4), round((ar_size(1)-96)/2)]);
asym_size = max(find(data(1:end-1,1,1)==0));
asym_percent = asym_size/size(data,1);
ar_size(1) = ar_size(1)-crop*2;
data = ifftshift(ifft(fftshift(data),[],1));
data = data(crop+1:end-crop,:,:);
data = ifftshift(fft(fftshift(data),[],1));

    %%
    for r_ind = 1:ar_size(end)
        k_space = (zeros(ar_size([1,2,3,4,5,6]),'single'));
        weights = (zeros(1,1,ar_size(3),ar_size(4),ar_size(5),ar_size(6),'single'));
        r1_rows = find(ar_ind(:,5)==r_ind);
%         ar_ind_r1 = ar_ind(r1_rows,:);
        ar_ind_r1 = ar_ind;
        
        replace = zeros(1,length(ar_ind_r1));
        counter = 0;
        counter2 = 0;
        % k_space(:,:,ar_ind(ind,1),ar_ind(ind,2),ar_ind(ind,4)) = data(:,:,ind);
        
        rWeights_tmp = rWeights(:,r_ind);
        for ind = 1:length(ar_ind_r1)
            
            if rWeights(ar_ind_r1(ind,6),r_ind) > weights(:,:,ar_ind_r1(ind,1),ar_ind_r1(ind,2),ar_ind_r1(ind,3),ar_ind_r1(ind,4))
                k_space(:,:,ar_ind_r1(ind,1),ar_ind_r1(ind,2),ar_ind_r1(ind,3),ar_ind_r1(ind,4)) = data(:,:,ar_ind_r1(ind,6));
                weights(:,:,ar_ind_r1(ind,1),ar_ind_r1(ind,2),ar_ind_r1(ind,3),ar_ind_r1(ind,4)) = rWeights_tmp(ar_ind_r1(ind,6));
%                 replace(ar_ind_r1(ind,6),r_ind) = 1;
%                 counter = counter +1;
%                 replaceW(counter,r_ind) = rWeights(ar_ind_r1(ind,6),r_ind);
            else
%                 counter2 = counter2 +1;
%                 skip(counter2,r_ind) = rWeights(ar_ind_r1(ind,6),r_ind);
            end
            
            
        end
    
        k_space = permute(k_space,[1,3,4,2,6,5]);
        kb = k_space(:,:,:,:,:,1);
        kx = k_space(:,:,:,:,:,2);
        ky = k_space(:,:,:,:,:,3);
        kz = k_space(:,:,:,:,:,4);
        clear k_space
        
        weights = repmat(weights, ar_size(1), ar_size(2), 1, 1, 1, 1);
        weights = permute(weights, [1,3,4,2,6,5]);
        weightsB = weights(:,:,:,:,:,1);
        weightsX = weights(:,:,:,:,:,2);
        weightsY = weights(:,:,:,:,:,3);
        weightsZ = weights(:,:,:,:,:,4);
        clear weights

        sampB = zeros(size(squeeze(kx(:,:,:,1,:))),'logical');
        sampX = zeros(size(squeeze(kx(:,:,:,1,:))),'logical');
        sampY = zeros(size(squeeze(kx(:,:,:,1,:))),'logical');
        sampZ = zeros(size(squeeze(kx(:,:,:,1,:))),'logical');

        sampB(abs(squeeze(kb(:,:,:,1,:)))>0) =1;
        sampX(abs(squeeze(kx(:,:,:,1,:)))>0) =1;
        sampY(abs(squeeze(ky(:,:,:,1,:)))>0) =1;
        sampZ(abs(squeeze(kz(:,:,:,1,:)))>0) =1;
    
        % Enforce Asymmetric Echo
        sampB(1:round(size(sampB,1)*asym_percent),:,:,:,:) = 0;
        sampX(1:round(size(sampX,1)*asym_percent),:,:,:,:) = 0;
        sampY(1:round(size(sampY,1)*asym_percent),:,:,:,:) = 0;
        sampZ(1:round(size(sampZ,1)*asym_percent),:,:,:,:) = 0;
        
        kb = bsxfun(@times,kb,permute(sampB,[1,2,3,5,4]));
        kx = bsxfun(@times,kx,permute(sampX,[1,2,3,5,4]));
        ky = bsxfun(@times,ky,permute(sampY,[1,2,3,5,4]));
        kz = bsxfun(@times,kz,permute(sampZ,[1,2,3,5,4]));
        
        weightsB = bsxfun(@times,weightsB,permute(sampB,[1,2,3,5,4]));
        weightsX = bsxfun(@times,weightsX,permute(sampX,[1,2,3,5,4]));
        weightsY = bsxfun(@times,weightsY,permute(sampY,[1,2,3,5,4]));
        weightsZ = bsxfun(@times,weightsZ,permute(sampZ,[1,2,3,5,4]));
        
        weightsB = squeeze(weightsB(:,:,:,1,:));
        weightsX = squeeze(weightsX(:,:,:,1,:));
        weightsY = squeeze(weightsY(:,:,:,1,:));
        weightsZ = squeeze(weightsZ(:,:,:,1,:));
        
        % ==========================================================================================
        % MAKE RESOLUTION ISOTROPIC (OR AS CLOSE AS POSSIBLE)
        % (This needs more work)
        % ==========================================================================================
        % Maximum allowable matrix size as contraint
        % ******************************************
        maxCols = 96;
        maxLins = 96;
        maxPars = 72;
        
        NImageCols = maxCols;
        NImageLins = maxLins;
        NImagePars = maxPars;
        
        % ******************************************
        % Current resolution
        ResCol = SG_Reveal4D_Data.param.ResCol;
        ResLin = SG_Reveal4D_Data.param.ResLin;
        ResPar = SG_Reveal4D_Data.param.ResPar;
        
        
%         % First make sure Col (FE) doesn't exceed maximum
%         NImageCols = min([SG_Reveal4D_Data.param.NImageCols, maxCols]);
%         padLength = (NImageCols - size(kb,1)) / 2;
%         if padLength >= 0
%             % pad
%             kb = padarray(kb,[padLength,0,0,0,0],0,'both');
%             kx = padarray(kx,[padLength,0,0,0,0],0,'both');
%             ky = padarray(ky,[padLength,0,0,0,0],0,'both');
%             kz = padarray(kz,[padLength,0,0,0,0],0,'both');
%             sampB = padarray(sampB,[padLength,0,0,0],0,'both');
%             sampX = padarray(sampX,[padLength,0,0,0],0,'both');
%             sampY = padarray(sampY,[padLength,0,0,0],0,'both');
%             sampZ = padarray(sampZ,[padLength,0,0,0],0,'both');
%             weightsB = padarray(weightsB,[padLength,0,0,0],0,'both');
%             weightsX = padarray(weightsX,[padLength,0,0,0],0,'both');
%             weightsY = padarray(weightsY,[padLength,0,0,0],0,'both');
%             weightsZ = padarray(weightsZ,[padLength,0,0,0],0,'both');
%         else
%             % crop
%             padLength = -padLength;
%             kb = kb(padLength+1:end-padLength,:,:,:,:);
%             kx = kx(padLength+1:end-padLength,:,:,:,:);
%             ky = ky(padLength+1:end-padLength,:,:,:,:);
%             kz = kz(padLength+1:end-padLength,:,:,:,:);
%             sampB = sampB(padLength+1:end-padLength,:,:,:);
%             sampX = sampX(padLength+1:end-padLength,:,:,:);
%             sampY = sampY(padLength+1:end-padLength,:,:,:);
%             sampZ = sampZ(padLength+1:end-padLength,:,:,:);
%             weightsB = weightsB(padLength+1:end-padLength,:,:,:);
%             weightsX = weightsX(padLength+1:end-padLength,:,:,:);
%             weightsY = weightsY(padLength+1:end-padLength,:,:,:);
%             weightsZ = weightsZ(padLength+1:end-padLength,:,:,:);
%         end
%         ResCol = ResCol/((size(kb,1)+2*padLength)/size(kb,1));
%         
%         
%         % Next, make in-plane resolution isotropic, Lin (PE) if possible
%         NImageLins = min([NImageCols, maxLins]);
%         padLength = (NImageLins - size(kb,2)) / 2;
%         if padLength >= 0
%             % pad
%             kb = padarray(kb,[0,padLength,0,0,0],0,'both');
%             kx = padarray(kx,[0,padLength,0,0,0],0,'both');
%             ky = padarray(ky,[0,padLength,0,0,0],0,'both');
%             kz = padarray(kz,[0,padLength,0,0,0],0,'both');
%             sampB = padarray(sampB,[0,padLength,0,0],0,'both');
%             sampX = padarray(sampX,[0,padLength,0,0],0,'both');
%             sampY = padarray(sampY,[0,padLength,0,0],0,'both');
%             sampZ = padarray(sampZ,[0,padLength,0,0],0,'both');
%             weightsB = padarray(weightsB,[0,padLength,0,0],0,'both');
%             weightsX = padarray(weightsX,[0,padLength,0,0],0,'both');
%             weightsY = padarray(weightsY,[0,padLength,0,0],0,'both');
%             weightsZ = padarray(weightsZ,[0,padLength,0,0],0,'both');
%         else
%             % crop
%             padLength = - padLength;
%             kb = kb(:,padLength+1:end-padLength,:,:,:);
%             kx = kx(:,padLength+1:end-padLength,:,:,:);
%             ky = ky(:,padLength+1:end-padLength,:,:,:);
%             kz = kz(:,padLength+1:end-padLength,:,:,:);
%             sampB = sampB(:,padLength+1:end-padLength,:,:);
%             sampX = sampX(:,padLength+1:end-padLength,:,:);
%             sampY = sampY(:,padLength+1:end-padLength,:,:);
%             sampZ = sampZ(:,padLength+1:end-padLength,:,:);
%             weightsB = weightsB(:,padLength+1:end-padLength,:,:);
%             weightsX = weightsX(:,padLength+1:end-padLength,:,:);
%             weightsY = weightsY(:,padLength+1:end-padLength,:,:);
%             weightsZ = weightsZ(:,padLength+1:end-padLength,:,:);
%         end
%         ResLin = ResLin/((size(kb,2)+2*padLength)/size(kb,2));     
%         
% %         % Finally, make Par (partitions/slices) isotropic if possible
% %         ResPar_tmp = min([ResCol, ResLin]);
% %         if (size(kb,3)*ResPar)/ResPar_tmp <= maxPars
% %             NImagePars = (size(kb,3)*ResPar)/ResPar_tmp;
% %             padLength = floor((NImagePars - size(kb,3)) / 2);
% %             if padLength >= 0
% %                 % pad
% %                 kb = padarray(kb,[0,0,padLength,0,0],0,'both');
% %                 kx = padarray(kx,[0,0,padLength,0,0],0,'both');
% %                 ky = padarray(ky,[0,0,padLength,0,0],0,'both');
% %                 kz = padarray(kz,[0,0,padLength,0,0],0,'both');
% %                 sampB = padarray(sampB,[0,0,padLength,0],0,'both');
% %                 sampX = padarray(sampX,[0,0,padLength,0],0,'both');
% %                 sampY = padarray(sampY,[0,0,padLength,0],0,'both');
% %                 sampZ = padarray(sampZ,[0,0,padLength,0],0,'both');
% %                 weightsB = padarray(weightsB,[0,0,padLength,0],0,'both');
% %                 weightsX = padarray(weightsX,[0,0,padLength,0],0,'both');
% %                 weightsY = padarray(weightsY,[0,0,padLength,0],0,'both');
% %                 weightsZ = padarray(weightsZ,[0,0,padLength,0],0,'both');
% %             else
% %                 % crop
% %                 kb = kb(:,:,padLength+1:end-padLength,:,:);
% %                 kx = kx(:,:,padLength+1:end-padLength,:,:);
% %                 ky = ky(:,:,padLength+1:end-padLength,:,:);
% %                 kz = kz(:,:,padLength+1:end-padLength,:,:);
% %                 sampB = sampB(:,:,padLength+1:end-padLength,:);
% %                 sampX = sampX(:,:,padLength+1:end-padLength,:);
% %                 sampY = sampY(:,:,padLength+1:end-padLength,:);
% %                 sampZ = sampZ(:,:,padLength+1:end-padLength,:);
% %                 weightsB = weightsB(:,:,padLength+1:end-padLength,:);
% %                 weightsX = weightsX(:,:,padLength+1:end-padLength,:);
% %                 weightsY = weightsY(:,:,padLength+1:end-padLength,:);
% %                 weightsZ = weightsZ(:,:,padLength+1:end-padLength,:);
% %             end
% %             ResPar = ResPar/((size(kb,3)+2*padLength)/size(kb,3));
% %             
% %         else
% %             ResPar_tmp = max([ResCol, ResLin]);
% %             if (size(kb,3)*ResPar)/ResPar_tmp <= maxPars
% %                 NImagePars = (size(kb,3)*ResPar)/ResPar_tmp;
% %                 padLength = floor((NImagePars - size(kb,3)) / 2);
% %                 if padLength >= 0
% %                     % pad
% %                     kb = padarray(kb,[0,0,padLength,0,0],0,'both');
% %                     kx = padarray(kx,[0,0,padLength,0,0],0,'both');
% %                     ky = padarray(ky,[0,0,padLength,0,0],0,'both');
% %                     kz = padarray(kz,[0,0,padLength,0,0],0,'both');
% %                     sampB = padarray(sampB,[0,0,padLength,0],0,'both');
% %                     sampX = padarray(sampX,[0,0,padLength,0],0,'both');
% %                     sampY = padarray(sampY,[0,0,padLength,0],0,'both');
% %                     sampZ = padarray(sampZ,[0,0,padLength,0],0,'both');
% %                     weightsB = padarray(weightsB,[0,0,padLength,0],0,'both');
% %                     weightsX = padarray(weightsX,[0,0,padLength,0],0,'both');
% %                     weightsY = padarray(weightsY,[0,0,padLength,0],0,'both');
% %                     weightsZ = padarray(weightsZ,[0,0,padLength,0],0,'both');
% %                 else
% %                     % crop
% %                     kb = kb(:,:,padLength+1:end-padLength,:,:);
% %                     kx = kx(:,:,padLength+1:end-padLength,:,:);
% %                     ky = ky(:,:,padLength+1:end-padLength,:,:);
% %                     kz = kz(:,:,padLength+1:end-padLength,:,:);
% %                     sampB = sampB(:,:,padLength+1:end-padLength,:);
% %                     sampX = sampX(:,:,padLength+1:end-padLength,:);
% %                     sampY = sampY(:,:,padLength+1:end-padLength,:);
% %                     sampZ = sampZ(:,:,padLength+1:end-padLength,:);
% %                     weightsB = weightsB(:,:,padLength+1:end-padLength,:);
% %                     weightsX = weightsX(:,:,padLength+1:end-padLength,:);
% %                     weightsY = weightsY(:,:,padLength+1:end-padLength,:);
% %                     weightsZ = weightsZ(:,:,padLength+1:end-padLength,:);
% %                 end
% %                 ResPar = ResPar/((size(kb,3)+2*padLength)/size(kb,3));
% %                 
% %                 
% %             else
%                 NImagePars = maxPars;
%                 padLength = floor((NImagePars - size(kb,3)) / 2);
%                 if padLength >= 0
%                     % pad
%                     kb = padarray(kb,[0,0,padLength,0,0],0,'both');
%                     kx = padarray(kx,[0,0,padLength,0,0],0,'both');
%                     ky = padarray(ky,[0,0,padLength,0,0],0,'both');
%                     kz = padarray(kz,[0,0,padLength,0,0],0,'both');
%                     sampB = padarray(sampB,[0,0,padLength,0],0,'both');
%                     sampX = padarray(sampX,[0,0,padLength,0],0,'both');
%                     sampY = padarray(sampY,[0,0,padLength,0],0,'both');
%                     sampZ = padarray(sampZ,[0,0,padLength,0],0,'both');
%                     weightsB = padarray(weightsB,[0,0,padLength,0],0,'both');
%                     weightsX = padarray(weightsX,[0,0,padLength,0],0,'both');
%                     weightsY = padarray(weightsY,[0,0,padLength,0],0,'both');
%                     weightsZ = padarray(weightsZ,[0,0,padLength,0],0,'both');
%                 else
%                     % crop
%                     kb = kb(:,:,padLength+1:end-padLength,:,:);
%                     kx = kx(:,:,padLength+1:end-padLength,:,:);
%                     ky = ky(:,:,padLength+1:end-padLength,:,:);
%                     kz = kz(:,:,padLength+1:end-padLength,:,:);
%                     sampB = sampB(:,:,padLength+1:end-padLength,:);
%                     sampX = sampX(:,:,padLength+1:end-padLength,:);
%                     sampY = sampY(:,:,padLength+1:end-padLength,:);
%                     sampZ = sampZ(:,:,padLength+1:end-padLength,:);
%                     weightsB = weightsB(:,:,padLength+1:end-padLength,:);
%                     weightsX = weightsX(:,:,padLength+1:end-padLength,:);
%                     weightsY = weightsY(:,:,padLength+1:end-padLength,:);
%                     weightsZ = weightsZ(:,:,padLength+1:end-padLength,:);
%                 end
%                 ResPar = ResPar/((size(kb,3)+2*padLength)/size(kb,3));
%                 
% %             end  
% %         end
%         
%         
%         
% %         % First make in-plane resolution isotropic:
% %         ResCol = SG_Reveal4D_Data.param.ResCol;
% %         ResLin = SG_Reveal4D_Data.param.ResLin;
% %         ResPar = SG_Reveal4D_Data.param.ResPar;
% %         
% %         NImageLins = min(SG_Reveal4D_Data.param.NImageLins, maxLin);
% %         padLength = (NImageLins - size(kb,2)) / 2;
% %         
% %         kb = padarray(kb,[0,padLength,0,0,0],0,'both');
% %         kx = padarray(kx,[0,padLength,0,0,0],0,'both');
% %         ky = padarray(ky,[0,padLength,0,0,0],0,'both');
% %         kz = padarray(kz,[0,padLength,0,0,0],0,'both');
% %         
% %         sampB = padarray(sampB,[0,padLength,0,0],0,'both');
% %         sampX = padarray(sampX,[0,padLength,0,0],0,'both');
% %         sampY = padarray(sampY,[0,padLength,0,0],0,'both');
% %         sampZ = padarray(sampZ,[0,padLength,0,0],0,'both');
% %         
% %         weightsB = padarray(weightsB,[0,padLength,0,0],0,'both');
% %         weightsX = padarray(weightsX,[0,padLength,0,0],0,'both');
% %         weightsY = padarray(weightsY,[0,padLength,0,0],0,'both');
% %         weightsZ = padarray(weightsZ,[0,padLength,0,0],0,'both');
% %         
% %         ResLin = ResCol;
% %         
% %         % Second, check if even padding can be acheived for partitions
% %         
% %         padLength1 = 0;
% %         padLength2 = 0;
% %         padLength3 = 0;
% %         
% %         if ResCol < ResPar
% %             tmp = round((ResPar/ResCol)*size(kb,3)) - size(kb,3);
% %             if mod(tmp,2) == 0
% %                 padLength3 = tmp/2;
% %             else
% %                 padLength3 = (tmp-1)/2;
% %             end
% %             % Update resolution
% %             ResPar = ResPar/((size(kb,3)+2*padLength3)/size(kb,3));
% %             % Pad arrays
% %             kb = padarray(kb,[0,0,padLength3,0,0],0,'both');
% %             kx = padarray(kx,[0,0,padLength3,0,0],0,'both');
% %             ky = padarray(ky,[0,0,padLength3,0,0],0,'both');
% %             kz = padarray(kz,[0,0,padLength3,0,0],0,'both');
% %             sampB = padarray(sampB,[0,0,padLength3,0],0,'both');
% %             sampX = padarray(sampX,[0,0,padLength3,0],0,'both');
% %             sampY = padarray(sampY,[0,0,padLength3,0],0,'both');
% %             sampZ = padarray(sampZ,[0,0,padLength3,0],0,'both');
% %             weightsB = padarray(weightsB,[0,0,padLength3,0],0,'both');
% %             weightsX = padarray(weightsX,[0,0,padLength3,0],0,'both');
% %             weightsY = padarray(weightsY,[0,0,padLength3,0],0,'both');
% %             weightsZ = padarray(weightsZ,[0,0,padLength3,0],0,'both');
% %             
% %         elseif ResCol > ResPar
% %             tmp1 = round((ResCol/ResPar)*size(kb,1)) - size(kb,1);
% %             tmp2 = round((ResLin/ResPar)*size(kb,2)) - size(kb,2);
% %             if mod(tmp1,2) == 0
% %                 padLength1 = tmp1/2;
% %                 if mod(tmp2,2) == 0
% %                     padLength2 = tmp2/2;
% %                 else
% %                     padLength2 = (tmp2-1)/2;
% %                 end
% %             else
% %                 padLength1 = (tmp1-1)/2;
% %                 if mod(tmp2,2) == 0
% %                     padLength2 = tmp2/2;
% %                 else
% %                     padLength2 = (tmp2-1)/2;
% %                 end
% %             end
% %             % Update resolution
% %             ResCol = ResCol/((size(kb,1)+2*padLength1)/size(kb,1));
% %             ResLin = ResLin/((size(kb,2)+2*padLength2)/size(kb,2));
% %             % Pad arrays
% %             kb = padarray(kb,[padLength1,padLength2,0,0,0],0,'both');
% %             kx = padarray(kx,[padLength1,padLength2,0,0,0],0,'both');
% %             ky = padarray(ky,[padLength1,padLength2,0,0,0],0,'both');
% %             kz = padarray(kz,[padLength1,padLength2,0,0,0],0,'both');
% %             sampB = padarray(sampB,[padLength1,padLength2,0,0],0,'both');
% %             sampX = padarray(sampX,[padLength1,padLength2,0,0],0,'both');
% %             sampY = padarray(sampY,[padLength1,padLength2,0,0],0,'both');
% %             sampZ = padarray(sampZ,[padLength1,padLength2,0,0],0,'both');
% %             weightsB = padarray(weightsB,[padLength1,padLength2,0,0],0,'both');
% %             weightsX = padarray(weightsX,[padLength1,padLength2,0,0],0,'both');
% %             weightsY = padarray(weightsY,[padLength1,padLength2,0,0],0,'both');
% %             weightsZ = padarray(weightsZ,[padLength1,padLength2,0,0],0,'both');
% %             
% %         end
            
        scanParam = SG_Reveal4D_Data.param;
        scanParam.NImageCols = size(kb,1);
        scanParam.NImageLins = size(kb,2);
        scanParam.NImagePars = size(kb,3);
        scanParam.ResCol = ResCol;
        scanParam.ResLin = ResLin;
        scanParam.ResPar = ResPar;

        
        
        dirName = [saveLocation,'\',name];
        if ~exist(dirName,'dir')
            mkdir(dirName);
        end
        
        % estimate acceleration rate
        % do not include asymmetric echo in rate calculation
        Rate_s = (size(sampB,2)*size(sampB,3)*size(sampB,4)) / sum(sum(sum(sampB(size(sampB,1)/2,:,:,:),2),3),4);
        Rate_w = (size(sampB,2)*size(sampB,3)*size(sampB,4)) / sum(sum(sum(weightsB(size(sampB,1)/2,:,:,:).^2,2),3),4);

%         Rate_s = prod(size(sampB)) / sum(sampB(:));
%         Rate_w = prod(size(sampB)) / sum(weightsB(:).^2);
        disp(['Acceleration rate...']);
        disp(['Samples: ',num2str(Rate_s)]);
        disp(['Weighted: ',num2str(Rate_w)]);
        
        save([dirName,'\',name,'_resp',num2str(r_ind),'_dataB'],'kb','-v7.3');
        save([dirName,'\',name,'_resp',num2str(r_ind),'_dataX'],'kx','-v7.3');
        save([dirName,'\',name,'_resp',num2str(r_ind),'_dataY'],'ky','-v7.3');
        save([dirName,'\',name,'_resp',num2str(r_ind),'_dataZ'],'kz','-v7.3');

        save([dirName,'\',name,'_resp',num2str(r_ind),'_sampB'],'sampB','-v7.3');
        save([dirName,'\',name,'_resp',num2str(r_ind),'_sampX'],'sampX','-v7.3');
        save([dirName,'\',name,'_resp',num2str(r_ind),'_sampY'],'sampY','-v7.3');
        save([dirName,'\',name,'_resp',num2str(r_ind),'_sampZ'],'sampZ','-v7.3');
        
        save([dirName,'\',name,'_resp',num2str(r_ind),'_weightsB'],'weightsB','-v7.3');
        save([dirName,'\',name,'_resp',num2str(r_ind),'_weightsX'],'weightsX','-v7.3');
        save([dirName,'\',name,'_resp',num2str(r_ind),'_weightsY'],'weightsY','-v7.3');
        save([dirName,'\',name,'_resp',num2str(r_ind),'_weightsZ'],'weightsZ','-v7.3');
        
        save([dirName,'\',name,'_resp',num2str(r_ind),'_scanParam'],'scanParam','-v7.3');
        save([dirName,'\',name,'_resp',num2str(r_ind),'_raw_header'],'raw_header','-v7.3');

    end
    
    clear kx ky kz
    
%     k_space = padarray(k_space,[padLength1,padLength+padLength2,padLength3,0,0,0],0,'both');
    k_space_n = kb;
    k_space_n(abs(k_space_n)==0) = NaN;
    k_space_n = mean(k_space_n,5,'omitnan');
%     k_space_n = mean(k_space_n,6,'omitnan');
    k_space_n(isnan(k_space_n))=0;
    k_space_n = squeeze(k_space_n);
    im = ifft3_shift(k_space_n);
    im = sos_combine(im);

    %%
    figure;
    imagesc(abs(im(:,:,round(size(im,3)/2)))); axis('off','image'); colormap('gray')

end



