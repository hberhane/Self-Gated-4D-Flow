classdef ReVEAL < hgsetget
% The ReVEAL Class reconstructs accelerated PCMRI data using the ReVEAL Method
%
% Methods:
%   ReVEAL: Constructor
%               Inputs: 
%
%   importData:Import the PC-MRI data
%               Inputs: fileName - The name of the file that contains the
%                        PC-MRI data.  The file should be "fileName.mat" in
%                        in the Matlab path
%
%   Method:
%               Inputs: 
%
%               Outputs:
%
%   Method:
%               Inputs: 
%
%               Outputs:
%
%   Method:
%               Inputs: 
%
%               Outputs:
%
% Properties
%   data
%
%   samplingPatterns
%
%   maxwellCorrections
%
%   options
%
%   parameters
%
%   scan_parameters
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich
% Written on:   2/28/2015
% Last update:  3/2/2015
%***************************************************************************
    
    % =====================================================================
    % Public Properties
    % =====================================================================
    properties
        data;
        samplingPatterns;
        maxwellCorrections;
        options;
        outputs;
        parameters;
        scan_parameters;
        saveDir;
    end
	
	% Private Properties
	properties(Access = protected, Hidden = true)
	end
    
    % =====================================================================
    % Public Methods
    % =====================================================================
    methods
        % =================================================================
        % Constructor Method
        function obj = ReVEAL()
            % Initialize ReVEALData object
            obj.data = ReVEALData();
            obj.options = ReVEALOptions();
            obj.outputs = ReVEALOutputs();
        end
        
        % =================================================================
        % Import Data
        importData(obj,fileName);
        
        % =================================================================
        % check data for errors
        dataChecks(obj);
        
        % =================================================================
        % Crop data along frequency encode dimension
        cropData(obj,crop_val);
        
        % =================================================================
        % Chunk data into smaller pieces to reconstruct separately
        chunkData(obj, dim, numChunks, overlap_prc);
        
        % =================================================================
        % Combine chunks into one image
        combineChunks(obj, dim);
        
        % =================================================================
        % Create dicoms from reconstructed images using dicom shell
        createDicomHeader(obj, dicomShellPath, saveName, medfilt, patientName);
       
        % =================================================================
        % Use Single Precision
        function useSingle(obj)
            obj.data.Yb = single(obj.data.Yb);
            obj.data.Yx = single(obj.data.Yx);
            obj.data.Yy = single(obj.data.Yy);
            obj.data.Yz = single(obj.data.Yz);
            obj.data.sampB = single(obj.data.sampB);
            obj.data.sampX = single(obj.data.sampX);
            obj.data.sampY = single(obj.data.sampY);
            obj.data.sampZ = single(obj.data.sampZ);
            obj.data.maxwellCorrX = single(obj.data.maxwellCorrX);
            obj.data.maxwellCorrY = single(obj.data.maxwellCorrY);
            obj.data.maxwellCorrZ = single(obj.data.maxwellCorrZ);
            obj.data.sensMaps = single(obj.data.sensMaps);
            obj.data.x0 = single(obj.data.x0);
%             obj.data.scanParam = single(obj.data.scanParam); % RA
            obj.options.ReVEALOpts.precision = 'single';
            
        end
        
        % =================================================================
        % Use GPU computation
        function useGPU(obj)
            if obj.options.isPlanar
                obj.data.Yb = gpuArray(obj.data.Yb);
                obj.data.Yx = gpuArray(obj.data.Yx);
                obj.data.Yy = gpuArray(obj.data.Yy);
                obj.data.Yz = gpuArray(obj.data.Yz);
                obj.data.sampB = gpuArray(obj.data.sampB);
                obj.data.sampX = gpuArray(obj.data.sampX);
                obj.data.sampY = gpuArray(obj.data.sampY);
                obj.data.sampZ = gpuArray(obj.data.sampZ);
                obj.data.maxwellCorrX = gpuArray(obj.data.maxwellCorrX);
                obj.data.maxwellCorrY = gpuArray(obj.data.maxwellCorrY);
                obj.data.maxwellCorrZ = gpuArray(obj.data.maxwellCorrZ);
                obj.data.sensMaps = gpuArray(obj.data.sensMaps);
%                 obj.data.scanParam = gpuArray(obj.data.scanParam); % RA
                obj.data.x0 = gpuArray(obj.data.x0);
            end
            obj.options.ReVEALOpts.compute = 'gpu';
        end
        
        % =================================================================
        % Pad kspace
        function padKspace(obj,padsize,dim)
%             if obj.options.isPlanar %RA
                if ~isempty(obj.data.sensMaps)
                    error('Pad kspace before sensitivity map estimation')
                else
                    if dim ==1
                        pad_array = [padsize];
                    elseif dim==2
                        pad_array = [0,padsize];
                    elseif dim==3
                        pad_array = [0,0,padsize];
                    end
                    obj.data.Yb = padarray(obj.data.Yb,pad_array,0,'both');
                    obj.data.Yx = padarray(obj.data.Yx,pad_array,0,'both');
                    obj.data.Yy = padarray(obj.data.Yy,pad_array,0,'both');
                    obj.data.Yz = padarray(obj.data.Yz,pad_array,0,'both');
                    obj.data.sampB = padarray(obj.data.sampB,pad_array,0,'both');
                    obj.data.sampX = padarray(obj.data.sampX,pad_array,0,'both');
                    obj.data.sampY = padarray(obj.data.sampY,pad_array,0,'both');
                    obj.data.sampZ = padarray(obj.data.sampZ,pad_array,0,'both');
                    %RA--added the following
                    obj.data.maxwellCorrX = padarray(obj.data.maxwellCorrX,pad_array,1,'both');
                    obj.data.maxwellCorrY = padarray(obj.data.maxwellCorrY,pad_array,1,'both');
                    obj.data.maxwellCorrZ = padarray(obj.data.maxwellCorrZ,pad_array,1,'both');
                end
%             end %RA
        end
        
        % =================================================================
        % Plot Sensitivity Maps
        function plotSense(obj)
            
            if ~obj.options.isPlanar
                figure
                clf
                title('Sensitivity Magnitudes')
                slice_center = round(size(obj.data.Yb,3)/2);
                for ind = 1:size(obj.data.Yb,4)
                    subplot(3,4,ind)
                    imagesc(squeeze(abs(obj.data.sensMaps(:,:,slice_center,ind,1))));
                    axis image
                    axis off
                end

                figure
                clf
                title('Sensitivity Magnitudes')
                for ind = 1:size(obj.data.Yb,4)
                    subplot(3,4,ind)
                    imagesc(squeeze(angle(obj.data.sensMaps(:,:,slice_center,ind,1))));
                    axis image
                    axis off
                end
            else
                figure
                clf
                title('Sensitivity Magnitudes')
                for ind = 1:size(obj.data.Yb,3)
                    subplot(3,4,ind)
                    imagesc(squeeze(abs(obj.data.sensMaps(:,:,ind,1))));
                    axis image
                    axis off
                end

                figure
                clf
                title('Sensitivity Magnitudes')
                for ind = 1:size(obj.data.Yb,3)
                    subplot(3,4,ind)
                    imagesc(squeeze(angle(obj.data.sensMaps(:,:,ind,1))));
                    axis image
                    axis off
                end
            end
            

            
        end
        
        % =================================================================
        % Estimate Sensitivy Maps
        function estimateSensMaps(obj,avg_all)
            
            if nargin < 2
                avg_all = 0;
            end
            % Check if the data is planar or volumetric
            if obj.options.isPlanar
                obj.estimateSensMaps2D();
            else
                obj.estimateSensMaps3D(avg_all);
            end            
        end
        
        % =================================================================
        % Reconstruct Data with ReVEAL
        function ReVEALRecon(obj,varargin)
            
            %Copy Options
            obj.options.GAMPOptB = obj.options.GAMPOpt.copy(obj.options.GAMPOpt);
            obj.options.GAMPOptX = obj.options.GAMPOpt.copy(obj.options.GAMPOpt);
            obj.options.GAMPOptY = obj.options.GAMPOpt.copy(obj.options.GAMPOpt);
            obj.options.GAMPOptZ = obj.options.GAMPOpt.copy(obj.options.GAMPOpt);
            
            if ~isempty(varargin)
                dir = varargin{1};
            else
                dir = 'x';
            end
            if obj.options.is1Dir
                obj.ReVEALRecon1Dir(dir);
            elseif ~obj.options.is1Dir && obj.options.isPlanar
                obj.ReVEALRecon3Dir();
            elseif ~obj.options.is1Dir && ~obj.options.isPlanar
%                 obj.ReVEALRecon4D_par();
                obj.ReVEALRecon4D();
            end
        end
             
        % =================================================================
        % Batch ReVEAL REcon
        function batchReVEALRecon(obj,fname,ID,joint)
            % Import data
            obj.importData(fname);

            % Do some data checks
            obj.dataChecks();
            
            % warning('cropping data')
            % obj.cropData(10);

            % Estimate the Sensitivity MAPs
            obj.estimateSensMaps2D();
            
            tmp_opts = obj.options.copy(obj.options);
            if ~joint
                % Reconstruct the data
                obj.ReVEALRecon('x');
                obj.options = tmp_opts.copy(obj.options);
                obj.ReVEALRecon('y');
                obj.options = tmp_opts.copy(obj.options);
                obj.ReVEALRecon('z');
            else
                obj.ReVEALRecon();
            end
            
            % Make a directory for the results
            curDir = pwd;
            path = getPath(fname);
            mkdir([path,'/',ID]);
            cd([path,'/',ID]);
            
            % Save The Results
            obj.saveData([fname(1:end-4),'_',ID]);
            obj.saveGIFs([fname(1:end-4),'_',ID]);
            cd(curDir);
        end
        
        % =================================================================
        % Estimate parameters after reconstruction
        function estimateParams(obj)
            
        end
        
        % =================================================================
        % Save the results
        function saveData( obj, fileName)
            
            % gather the outputs before saving
            obj.outputs.xHatb = gather(obj.outputs.xHatb);
            obj.outputs.xHatx = gather(obj.outputs.xHatx);
            obj.outputs.xHaty = gather(obj.outputs.xHaty);
            obj.outputs.xHatz = gather(obj.outputs.xHatz);
            obj.outputs.thetaX = gather(obj.outputs.thetaX);
            obj.outputs.thetaY = gather(obj.outputs.thetaY);
            obj.outputs.thetaZ = gather(obj.outputs.thetaZ);
            obj.outputs.thetaMAPX = gather(obj.outputs.thetaMAPX);
            obj.outputs.thetaMAPY = gather(obj.outputs.thetaMAPY);
            obj.outputs.thetaMAPZ = gather(obj.outputs.thetaMAPZ);
            obj.outputs.vX = gather(obj.outputs.vX);
            obj.outputs.vY = gather(obj.outputs.vY);
            obj.outputs.vZ = gather(obj.outputs.vZ);
            
            % Copy outputs and options for saving
            outputs = ReVEALOutputs();
            outputs = obj.outputs;

            options = ReVEALOptions();
            options = obj.options;
            
            % Get rid variables taking up memory
            options.GAMPOptB.xhat0 = [];
            options.GAMPOptB.xvar0 = [];
            options.GAMPOptB.shat0 = [];
            options.GAMPOptB.svar0 = [];
            options.GAMPOptB.xhatPrev0 = [];
            
            options.GAMPOptX.xhat0 = [];
            options.GAMPOptX.xvar0 = [];
            options.GAMPOptX.shat0 = [];
            options.GAMPOptX.svar0 = [];
            options.GAMPOptX.xhatPrev0 = [];
            
            options.GAMPOptY.xhat0 = [];
            options.GAMPOptY.xvar0 = [];
            options.GAMPOptY.shat0 = [];
            options.GAMPOptY.svar0 = [];
            options.GAMPOptY.xhatPrev0 = [];
            
            options.GAMPOptZ.xhat0 = [];
            options.GAMPOptZ.xvar0 = [];
            options.GAMPOptZ.shat0 = [];
            options.GAMPOptZ.svar0 = [];
            options.GAMPOptZ.xhatPrev0 = [];
            
            % save outputs and optio
            s = warning('error', 'MATLAB:save:sizeTooBigForMATFile');
            warning('error', 'MATLAB:save:sizeTooBigForMATFile');
            % Normal save
            try
                save([fileName,'_results'],'outputs');
                save([fileName,'_options'],'options');
            % If too big, save v7.3
            catch
                save([fileName,'_results'],'outputs', '-v7.3');
                save([fileName,'_options'],'options');
            end
            warning(s);
                
        end
        
        % =================================================================
        % Save GIFs of the results
        function saveGIFs( obj, fileName,varargin)
            if isempty(varargin)
                fps = 10;
                clip = 2;
            else
                fps = varargin{1};
                clip = varargin{2};
            end
            create_GIF(gather(abs(obj.outputs.xHatb)),[fileName,'_pB'],fps,clip);
            if ~isempty(obj.outputs.xHatx)
                create_GIF(gather(abs(obj.outputs.xHatx)),[fileName,'_pX'],fps,clip);
                create_GIF(gather(norm_velocity(obj.outputs.thetaX)),[fileName,'_thetaX'],fps,1);
                create_GIF(gather(obj.outputs.vX),[fileName,'_vX'],fps,1);
            end
            if ~isempty(obj.outputs.xHaty)
                create_GIF(gather(abs(obj.outputs.xHaty)),[fileName,'_pY'],fps,clip);
                create_GIF(gather(norm_velocity(obj.outputs.thetaY)),[fileName,'_thetaY'],fps,1);
                create_GIF(gather(obj.outputs.vY),[fileName,'_vY'],fps,1);
            end
            if ~isempty(obj.outputs.xHatz)
                create_GIF(gather(abs(obj.outputs.xHatz)),[fileName,'_pZ'],fps,clip);
                create_GIF(gather(norm_velocity(obj.outputs.thetaZ)),[fileName,'_thetaZ'],fps,1);
                create_GIF(gather(obj.outputs.vZ),[fileName,'_vZ'],fps,1);
            end
        end
        
        function permuteImage( obj, pattern )
            obj.outputs.xHatb = permute(obj.outputs.xHatb, pattern);
            obj.outputs.xHatx = permute(obj.outputs.xHatx, pattern);
            obj.outputs.xHaty = permute(obj.outputs.xHaty, pattern);
            obj.outputs.xHatz = permute(obj.outputs.xHatz, pattern);
            obj.outputs.thetaX = permute(obj.outputs.thetaX, pattern);
            obj.outputs.thetaY = permute(obj.outputs.thetaY, pattern);
            obj.outputs.thetaZ = permute(obj.outputs.thetaZ, pattern);
            obj.outputs.vX = permute(obj.outputs.vX, pattern);
            obj.outputs.vY = permute(obj.outputs.vY, pattern);
            obj.outputs.vZ = permute(obj.outputs.vZ, pattern); 
        end
        
        ReVEALRecon4D_1V(obj)
        
    end    
    
    % =====================================================================
    % Private Methods
    % =====================================================================
    methods(Access = protected, Hidden = true)
        
        % Method to Estimate 2D sensitivity Maps
        estimatSensMaps2D(obj);
        
        % Method to Estimate 3D sensitivity Maps
        estimatSensMaps3D(obj);
        
        % Reconstruct Data with ReVEAL for 1 velocity directions
        ReVEALRecon1Dir(obj,dir);
        
        % Reconstruct Data with ReVEAL for 3 velocity directions
        ReVEALRecon3Dir(obj);
        
        % Reconstruct Data with ReVEAL for 3 velocity directions
        ReVEALRecon4D(obj);
        
    end
end

