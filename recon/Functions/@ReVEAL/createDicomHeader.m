function [] = createDicomHeader(Obj, dicomShellPath, saveName, medfilt, patientName)

% addpath(genpath('C:\Users\tesla\Downloads\aaron_resp_data'))
% addpath(genpath('C:\Users\tesla\Documents\MATLAB\Adam'))

xHatb = Obj.outputs.xHatb;
xHatx = Obj.outputs.xHatx;
xHaty = Obj.outputs.xHaty;
xHatz = Obj.outputs.xHatz;

% **************************************************************************************************
% % TEST to denoise complex difference
% % Seems promising, experiement later with block-hankel structure
% c_x = xHatb-xHatx;
% c_y = xHatb-xHaty;
% c_z = xHatb-xHatz;
% c_b = xHatb+xHatx;
% 
% A = cat(4, c_x, c_y, c_z, c_b);
% A = reshape(A, [size(c_x,1)*size(c_x,2)*size(c_x,3), size(A,4)]);
% 
% [V,D] = eig(A'*A);
% % Sort eigenvalues and eigenvectors in descending order
% [~,ind] = sort(diag(D),'descend');
% S = diag(D); S = S(ind);
% % Sorted eigenvectors
% Vs = V(:,ind);
% % Singular vectors
% US  = A*Vs;
% US(:,11:end) = 0;
% 
% A_new = US*Vs';
% 
% I = reshape(A_new, [size(c_x,1),size(c_x,2),size(c_x,3),size(c_x,4),4]);
% 
% I_x = I(:,:,:,:,1);
% I_y = I(:,:,:,:,2);
% I_z = I(:,:,:,:,3);
% I_b = I(:,:,:,:,4);
% 
% b = (I_b + I_x)/2;
% x = -I_x + b;
% y = -I_y + b;
% z = -I_z + b;
% 
% blah = b - xHatb;
% 
% theta_x = angle(xHatx.*conj(xHatb));
% theta_x_new = angle(x.*conj(b));
% 
% 
% figure;
% for j = 1 : 10
% for i = 1 : 20
%     subplot(121); imagesc(squeeze(abs(xHatb(:,:,45,i).^0.8))); axis('off','image'); colormap('gray');
%     subplot(122); imagesc(squeeze(abs(b(:,:,45,i).^0.8))); axis('off','image'); colormap('gray'); pause(0.1);
% end
% end
% 
% figure; imagesc(squeeze(abs(b(:,:,45,5).^0.8))); axis('off','image'); colormap('gray');
% 
% figure; imagesc(squeeze((theta_x(:,:,45,5))),[-pi,pi]); axis('off','image'); colormap('gray');
% figure; imagesc(squeeze((theta_x_new(:,:,45,5))),[-pi,pi]); axis('off','image'); colormap('gray');
% 
% figure;
% for i = 1 : 20
% imagesc(squeeze((theta_x_new(:,:,45,i))),[-pi,pi]); axis('off','image'); colormap('gray'); pause(0.2); 
% end
% 
% % **************************************************************************************************

% cDif = sqrt((xHatb - xHatx).^2 + (xHatb - xHaty).^2 + (xHatb - xHatz).^2);
% cDif = abs(cDif);
% Dilate using gaussian kernal
k_1(:,1,1) = [1,2,1];
k_2(1,:,1) = [1,2,1];
k_3(1,1,:) = [1,2,1];
kernal = bsxfun(@times, bsxfun(@times, k_1, k_2), k_3);
kernal = kernal/sum(kernal(:));

% for i = 1 : size(cDif, 4)
%     cDif(:,:,:,i) = convn(cDif(:,:,:,i),kernal,'same');
% end
% 
% cDif = max(bsxfun(@times, max(cDif,[],4), cDif),[],4);
% 
% cDif = max(abs(cDif),[], 4);
% Mag_tilda = cDif;
% cDif = mean(abs(cDif),4);
% cDif = repmat(cDif, [1,1,1,20]);
% Mag_tilda = abs(xHatb).*abs(cDif);
% Mag_tilda = mean(Mag_tilda, 4);
% Mag_tilda = repmat(Mag_tilda, [1,1,1,20]);

im4D.M      = abs(Obj.outputs.xHatb);
% im4D.C      = Mag_tilda;
im4D.X      = Obj.outputs.thetaX;
im4D.Y      = Obj.outputs.thetaY;
im4D.Z      = Obj.outputs.thetaZ;
scanParam   = Obj.outputs.scanParam;

% Only for SGvPT
im4D.M = padarray(im4D.M, [20, 0, 0, 0], 'pre');
% im4D.C = padarray(im4D.C, [20, 0, 0, 0], 'pre');
im4D.X = padarray(im4D.X, [20, 0, 0, 0], 'pre');
im4D.Y = padarray(im4D.Y, [20, 0, 0, 0], 'pre');
im4D.Z = padarray(im4D.Z, [20, 0, 0, 0], 'pre');

im4D.M = flip(flip(im4D.M, 2), 3);
% im4D.C = flip(flip(im4D.C, 2), 3);
im4D.X = flip(flip(im4D.X, 2), 3);
im4D.Y = flip(flip(im4D.Y, 2), 3);
im4D.Z = flip(flip(im4D.Z, 2), 3);

% Crop oversampling if necessary
crop1 = round((size(im4D.M,1) - scanParam.NImageCols)/2);
crop2 = round((size(im4D.M,2) - scanParam.NImageLins)/2);
crop3 = round((size(im4D.M,3) - scanParam.slThck/scanParam.ResPar)/2);

% crop1 = (96 - 96)/2;
% crop2 = (96 - 96)/2;
% crop3 = (72 - 56)/2;

im4D.M = im4D.M(crop1+1:end-crop1,crop2+1:end-crop2,crop3+1:end-crop3,:);
% im4D.C = im4D.C(crop1+1:end-crop1,crop2+1:end-crop2,crop3+1:end-crop3,:);
im4D.X = im4D.X(crop1+1:end-crop1,crop2+1:end-crop2,crop3+1:end-crop3,:);
im4D.Y = im4D.Y(crop1+1:end-crop1,crop2+1:end-crop2,crop3+1:end-crop3,:);
im4D.Z = im4D.Z(crop1+1:end-crop1,crop2+1:end-crop2,crop3+1:end-crop3,:);

% make Y velocity component negative
% why is this necessary? Because stupid reasons...
im4D.Y = -im4D.Y;
im4D.Z = -im4D.Z;



% % Background phase correction using WRLS+ARTO
opt.fit = 'WRLS+ARTO';
opt.pOrd = 2;
opt.lam = 50;
opt.mTh = 0.04;
opt.midFOVFrac = 0.5;
opt.selectInit = 0;
opt.Kmax = 2;
opt.tau = 3;
opt.delta = 2;
opt.gmmComp = 3;
for i = 1 : 3
    switch i
        case 1
            phase_cine = im4D.X;
        case 2
            phase_cine = im4D.Y;
        case 3
            phase_cine = im4D.Z;
    end
    
    % Express phase as a percentage of the range
    phase_cine = (double(phase_cine)/pi);
    
    % Pre-processing before background phase correction
    % Time-averaged phase
    phi = mean(phase_cine,4);
    % Temporal standard deviation
    sigma = std(phase_cine, 0, 4);
    % Magnitude thresholding
    mag_avg = mean(im4D.M, 4);
    M = ones(size(mag_avg));
    M(mag_avg <= opt.mTh*max(mag_avg(:))) = 0;
    M = logical(M);
    
    % Background phase correction using WRLS+ARTO
    [phi_Cor, cMap] = wrls_arto(phi, sigma, M, opt);
    figure; imagesc(squeeze(cMap(:,:,28)),[-0.05 0.05]); axis('off','image'); colormap('jet')
    figure; imagesc(squeeze(phi(:,:,28)),[-0.05 0.05]); axis('off','image'); colormap('jet')
    figure; imagesc(squeeze(phi_Cor(:,:,28)),[-0.05 0.05]); axis('off','image'); colormap('jet')
    
    phase_cine_corrected = bsxfun(@minus, phase_cine, cMap);
    phase_cine_corrected = single(phase_cine_corrected*pi);
    
    switch i
        case 1
            im4D.X = phase_cine_corrected;
        case 2
            im4D.Y = phase_cine_corrected;
        case 3
            im4D.Z = phase_cine_corrected;
    end
end
    

% ------------------------------------------------------------------------
if medfilt
    % Apply median filter (optional) to improve peak velocity quantification
    for encoding = 1 : 3
        switch encoding
            case 1
                im = im4D.X;
            case 2
                im = im4D.Y;
            case 3
                im = im4D.Z;
        end
        tmp = zeros(size(im4D.X));
        
        for i = 2 : size(im4D.X,1)-1
            for j = 2 : size(im4D.X,2)-1
                for k = 2 : size(im4D.X,3)-1
                    
                    for c = 1 : size(im4D.X,4)
                        region(1) = im(i,j,k,c);
                        region(2) = im(i-1,j,k,c);
                        region(3) = im(i+1,j,k,c);
                        region(4) = im(i,j-1,k,c);
                        region(5) = im(i,j+1,k,c);
                        region(6) = im(i,j,k-1,c);
                        region(7) = im(i,j,k+1,c);
                        
                        tmp(i,j,k,c) = median(region);
                    end
                end
            end
        end
        
        tmp(1,:,:) = im(1,:,:);
        tmp(end,:,:) = im(end,:,:);
        tmp(:,1,:) = im(:,1,:);
        tmp(:,end,:) = im(:,end,:);
        tmp(:,:,1) = im(:,:,1);
        tmp(:,:,end) = im(:,:,end);
        
        switch encoding
            case 1
                im4D.X = tmp;
            case 2
                im4D.Y = tmp;
            case 3
                im4D.Z = tmp;
        end
    end
end
% -------------------------------------------------------------------------


% *************************************************************************
% Heart Segmentation
% mask = logical(zeros(size(im4DAvg)));
% ind = flip(1:size(im4DAvg,3));
% for i = 1 : size(im4DAvg,3)
%     i
%     mask(:,:,i) = roipoly(im4DAvg(:,:,i));
% end
% 
% % mask = permute(mask,[3,2,1]);
% 
% im4D.M = bsxfun(@times, im4D.M, mask);
% im4D.X = bsxfun(@times, im4D.X, mask);
% im4D.Y = bsxfun(@times, im4D.Y, mask);
% im4D.Z = bsxfun(@times, im4D.Z, mask);
% *************************************************************************

% *************************************************************************
% Interpolate images
% im4DMInt = zeros(2*size(im4D.M,1),2*size(im4D.M,2),2*size(im4D.M,3),size(im4D.M,4));
% im4DXInt = zeros(2*size(im4D.M,1),2*size(im4D.M,2),2*size(im4D.M,3),size(im4D.M,4));
% im4DYInt = zeros(2*size(im4D.M,1),2*size(im4D.M,2),2*size(im4D.M,3),size(im4D.M,4));
% im4DZInt = zeros(2*size(im4D.M,1),2*size(im4D.M,2),2*size(im4D.M,3),size(im4D.M,4));
% 
% for i = 1 : 20
%     im4DMInt(:,:,:,i) = imresize3(im4D.M(:,:,:,i),2);
%     im4DXInt(:,:,:,i) = imresize3(im4D.X(:,:,:,i),2);
%     im4DYInt(:,:,:,i) = imresize3(im4D.Y(:,:,:,i),2);
%     im4DZInt(:,:,:,i) = imresize3(im4D.Z(:,:,:,i),2);
% end
% 
% im4D.M = im4DMInt;
% im4D.X = im4DXInt;
% im4D.Y = im4DYInt;
% im4D.Z = im4DZInt;
% 
% clear im4DMInt im4DXInt im4DYInt im4DZInt
% *************************************************************************

% Make new "magnitude" images my multiplying Mag w/ complex difference
% Idea is to make segmentation easier in software
% im4D.M = im4D.C;

% Scale to 12-bit (0 to 4095)
% im4D.M = abs(im4D.M .^ 0.8); 
im4D.M = im4D.M - min(im4D.M(:)); im4D.M = im4D.M / max(im4D.M(:)); im4D.M = im4D.M * 4095;
im4D.X = (im4D.X+pi)/(2*pi); im4D.X = im4D.X * 4095;
im4D.Y = (im4D.Y+pi)/(2*pi); im4D.Y = im4D.Y * 4095;
im4D.Z = (im4D.Z+pi)/(2*pi); im4D.Z = im4D.Z * 4095;

% Load dicom shells
shellDir = dir(dicomShellPath);
fname{1} = [dicomShellPath,shellDir(3).name,'\'];
fname{2} = [dicomShellPath,shellDir(4).name,'\'];
fname{3} = [dicomShellPath,shellDir(5).name,'\'];
fname{4} = [dicomShellPath,shellDir(6).name,'\'];

% % Load segmentation mask (optional)
% mask = load('C:\Users\spc\Documents\MATLAB\wholeheartMask_Tot_Thresh.mat'); mask = mask.mask_tot_thresh;
% mask = flip(flip(mask, 2), 3);
% mask = mask(crop1+1:end-crop1,crop2+1:end-crop2,crop3+1:end-crop3,:);

for ind = 1 : 4
[a_0, t_info, N_S, fnames] = Read_Dicom_Folder_v3(fname{ind});

    switch ind
        case 1
            tmp = uint16(im4D.M);
        case 2
            tmp = uint16(im4D.X);
        case 3
            tmp = uint16(im4D.Y);
        case 4
            tmp = uint16(im4D.Z);
    end
    
%     tmp = tmp .* repmat(uint16(mask), [1,1,1,20]);
    
    nCPhases = 20; % Need to add this as field to "scanParam"
    seriesNumber = t_info(1,1).SeriesNumber;
    NominalInterval = 1000*scanParam.meanRR;
    TriggerTime = [0:1:nCPhases-1]*(NominalInterval/nCPhases);
    ResCol = scanParam.ResCol;
    ResLin = scanParam.ResLin;
    ResPar = scanParam.ResPar;
%     FOV = [
    
    
    if ~exist([saveName,'\',num2str(seriesNumber)],'dir')
        mkdir([saveName,'\',num2str(seriesNumber)]);
    end

    cd([saveName,'\',num2str(seriesNumber)]);
    instance = 0;
    for j = 1 : N_S
        for i = 1 : nCPhases
            instance = instance+1;
            fileName = ['MR00',num2str(instance),'.dcm'];
        
            header = t_info(j,1);
            header.NominalInterval = NominalInterval;
            header.TriggerTime = TriggerTime(i);
            header.CardiacNumberOfImages = nCPhases;
            header.InstanceNumber = instance;
            header.SeriesNumber = seriesNumber;
            header.ImageComments = 'SGvPT';
%             header.PixelSpacing = [ResCol; ResLin];
%             header.SliceThickness = ResPar;
%             header.Private_0051_100c = 
            
            header.PatientName.FamilyName = patientName;
        
            currentIm = squeeze(tmp(:,:,j,i));
            dicomwrite(currentIm, fileName, header, 'CreateMode', 'copy','WritePrivate',true);
        end
    end
    
end

end












