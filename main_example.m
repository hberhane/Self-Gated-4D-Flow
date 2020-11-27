%==================================================================================================%
%\\\____________________________________________________________________________________________///%
%[  Main script
%[  Run Self-gated 4D flow binning with ReVEAL4D reconstruction
%[
%[________________________________________________________________________________________________]%
%///--------------------------------------------------------------------------------------------\\\%
clear;  close all;
addpath(genpath('..\Self-Gated-4D-Flow'));
%==================================================================================================%


% Options
%---------------------------------------------------------------------------------------------------
opt = init();
%************
opt.sigExt      = 'PCA';        % Signal extraction method, 'PCA' , 'SSAFARI', 'PT', or 'PTSola'
opt.rSoft       = 1;            % Soft-gating/weighting or not on respiratory dim
opt.nRPhases    = 1;            % # of respiratory phases to bin
opt.nCPhases    = 20;           % base # of cardiac phases to bin
opt.viewSharing = 0;            % view sharing to increase # of cardiac phases (x2)
opt.recStart    = 'start';     % 'start', 'center', 'end': where to begin trimming for "recTime"
opt.expOrd      = 4;
opt.respEff     = 50;
opt.flow        = 1;            % Flow dataset or not
opt.rampSS      = 50;           % Trim # self-gating samples during steady state ramp
opt.recTime     = 300;
opt.sgI         = 9;


% Load raw data (.dat) file
%---------------------------------------------------------------------------------------------------
opt.fName = which('meas_MID00396_FID18723_4D_FLOW_WIP_retro_ePAT_aorta.dat');
[Dat, ~, param, ~, rawData] = readWrapper4(opt);
rawData = rawData{2};
param.TR = param.TR * 1000;


% Bin the data
%---------------------------------------------------------------------------------------------------
FlowSG4D_Outputs_SG = binSelfGating(Dat, rawData, param, opt);


% Construct array
%---------------------------------------------------------------------------------------------------
name = 'example_binned';
saveLocation = '.\example_dataset';
raw_header = rawData.hdr;
[resCol,resLin,resPar,rate_samps,rate_weights,NImageCols,NImageLins,NImagePars] = constructArrayWrapper(FlowSG4D_Outputs_SG, opt, name, saveLocation, param, raw_header);
addpath(genpath(saveLocation));
       

% ReVEAL4D reconstruction
%---------------------------------------------------------------------------------------------------
% Find files
Fname = [name, '_resp1_'];
saveDir = fullfile(pwd,'recon');

% Initialize class and import data
op = ReVEAL();
op.saveDir = saveDir;

% Set some options
op.options.ReVEALOpts.gamma = 0.95;
op.options.ReVEALOpts.lambda0 = 1.5;
op.options.ReVEALOpts.nit = 8;
op.options.ReVEALOpts.uniform_var = 1;
op.options.ReVEALOpts.L1 = 1;
op.options.ReVEALOpts.MAP = 1;
op.options.GAMPOpt.nit = 40;
op.options.GAMPOpt.verbose = 1;
op.options.GAMPOpt.tol = 1e-4;
op.options.ReVEALOpts.wvar = 1e-10;

% Import data
op.importData([Fname,'data.mat']);
op.options.ReVEALOpts.sigmaSq = 0.01;

% Find sensitivity maps and set some other things
op.useSingle();
op.estimateSensMaps();
% gpuDevice(1); % GPU(1) or GPU(2)
% op.useGPU();
% op.cropData(5);

if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% Recon.
op.ReVEALRecon();


% Saving data and gifs
%---------------------------------------------------------------------------------------------------
op_acq = op;
op_interp = op;
% Acquired resolution
% *******************
fps = 10;
clip = 1.5;
saveDir_acq = [saveDir, '\Acquired'];
saveName = [fname{ind},'ReVEAL4D_Acq_','lam_',num2str(op.options.ReVEALOpts.lambda0),'_sSq_',num2str(op.options.ReVEALOpts.sigmaSq),'_gamma_',num2str(op.options.ReVEALOpts.gamma)];
saveName = strrep(saveName, '.', 'p');
if ~exist(saveDir_acq, 'dir')
    mkdir(saveDir_acq);
end
op_acq.saveData([saveDir_acq,'\',saveName])
save_gif_wrapper(op_acq, saveDir_acq, saveName, fps, clip);

%===================================================================================================
%===================================================================================================

 