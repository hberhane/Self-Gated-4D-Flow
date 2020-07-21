clear
close all;

addpath(genpath('D:\MRI Data\Aaron\4DSG MRM'));                                   % data
addpath(genpath('D:\MRI Data\Aaron\Joint Cine Flow (JCF)'));                      % data
addpath(genpath('D:\MRI Data\Aaron\CardioStepper4D'));                      % data
addpath(genpath('D:\MRI Data\Aaron\4DFlow SG vs PT'));                      % data
addpath(genpath('C:\Users\spc\Documents\MATLAB\Self_gating'));                    % self-gating code
addpath(genpath('C:\Users\spc\Documents\MATLAB\RACommons\readScannerData'));      % readwrapper
addpath(genpath('C:\Users\spc\Documents\MATLAB\RACommons\misc'));                 % misc
addpath(genpath('C:\Users\spc\Documents\MATLAB\Adam'));                           % misc
addpath(genpath('D:\MRI Data\Aaron\NCH FeraHeme 4D Flow (NCHFH4D)'));
addpath(genpath('D:\MRI Data\Aaron\Sola'));
addpath(genpath('D:\MRI Data\Aaron\TMP'));


%% =================================================================================================
% Data load
% ==================================================================================================
X = {    
        'meas_MID00049_FID18004_BEAT_FQ_SG_785A_WholeHeart_20190314.dat',                           'YL',           'VE',  19;  % 1
        'meas_MID00036_FID01782_BEAT_FQ_SG_785A_WholeHeart_20190314.dat',                           'CC',           'VE',  35;  % 2
        'meas_MID00041_FID03981_BEAT_FQ_SG_785A_WholeHeart_20190314.dat',                           'WG',           'VE',  23;  % 3
        'meas_MID00308_FID04489_BEAT_FQ_SG_785A_WholeHeart_20190314.dat',                           'NP',           'VE',  29;  % 4
        'meas_MID00039_FID08581_BEAT_FQ_SG_785A_WholeHeart_20190314.dat',                           'JS',           'VE',  15;  % 5
        'meas_MID00265_FID10798_BEAT_FQ_SG_785A_WholeHeart_20190314.dat',                           'CCT',          'VE',  45;  % 6
        'meas_MID00102_FID18057_BEAT_FQ_SG_785A_WholeHeart_20190314.dat',                           'AP',           'VE',  55;  % 7
        'meas_MID00094_FID18171_BEAT_FQ_SG_785A_WholeHeart_20190314.dat',                           'SL',           'VE',  24;  % 8
        'meas_MID00164_FID18241_BEAT_FQ_SG_785A_WholeHeart_20190314.dat',                           'YP',           'VE',  23;  % 9
        'meas_MID00436_FID25192_BEAT_FQ_SG_785A_WholeHeart_20190314_FF.dat',                        'FF',           'VE',  100; % 10   
        'meas_MID00307_FID21660_BEAT_FQ_SG_785A_WholeHeart_20190314.dat',                           'YL_R1',        'VE',  100; % 11
        'meas_MID00312_FID21665_BEAT_FQ_SG_785A_WholeHeart_20190314_Repeat2.dat',                   'YL_R2',        'VE',  100; % 12
        'meas_MID00317_FID21670_BEAT_FQ_SG_785A_WholeHeart_20190314_Repeat3.dat',                   'YL_R3',        'VE',  100; % 13
        'meas_MID00322_FID21675_BEAT_FQ_SG_785A_WholeHeart_20190314_Repeat4.dat',                   'YL_R4',        'VE',  100; % 14
        'meas_MID00327_FID21680_BEAT_FQ_SG_785A_WholeHeart_20190314_Repeat5.dat',                   'YL_R5',        'VE',  100; % 15
        'meas_MID00400_FID63895_BEAT_FQ_sg_Aorta_1min_FA_7.dat',                                    'BH_Ao',        'VE',  100; % 16
        'meas_MID00393_FID63888_BEAT_FQ_sg_whole_heart_5min_FA_5.dat',                              'FA_5',         'VE',  100; % 17
        'meas_MID00392_FID63887_BEAT_FQ_sg_whole_heart_5min_FA_7.dat',                              'FA_7',         'VE',  100; % 18
        'meas_MID00340_FID04812_CV_3dSSFP_SG_20190821_2min.dat',                                    'Cine_Phantom',   'VE',  100; % 19
        'meas_MID00413_FID04885_CV_3dSSFP_SG_20190821_2min.dat',                                    'Cine_RA',      'VE',  100; % 20
        'meas_MID00319_FID63814_OSU_NOISE_4CH.dat',                                                 'Cine_RA',      'VE',  100; % 21
        'meas_MID00078_FID67757_BEAT_FQ_sg_whole_heart_8min.dat',                                   'NCH_FeraHeme',   'VE',  100; % 22
        'meas_MID00042_FID73759_CV_ssfp3d_sg_20190904_FA90.dat',                                    'Cine_Avanto_FA90',   'VD',  100; % 23
        'meas_MID00059_FID74750_CV_ssfp3d_sg_20190904_FA90_LowSAR_Pulse.dat',                       'Cine_RA',      'VD',  100; % 24
        'meas_MID00068_FID74759_CV_ssfp3d_sg_20190904_FA90_LowSAR_Pulse_exercise.dat',              'Cine_RA_Exc',   'VD',  100; % 25
        'meas_MID00056_FID74747_BEAT_FQ_SG_whole_heart_pre.dat',                                    'Flow_RA',          'VD',  100; % 26
        'meas_MID00067_FID74758_BEAT_FQ_SG_whole_heart_exercise.dat',                               'Flow_RA_Exc',           'VD',  100; % 27
        'meas_MID00389_FID21409_BEAT_FQ_SG_whole_heart.dat',                                        'PT_KS',           'VD',  100; % 28
        'meas_MID00464_FID23436_BEAT_FQ_SG_whole_heart.dat',                                        'PT_NP',           'VD',  100; % 29
        'meas_MID00082_FID94301_BEAT_FQ_sg_whole_heart_8min.dat',                                   'NCH_FeraHeme_Mky',     'VE',  100; % 30
        'meas_MID00128_FID02389_BEAT_FQ_20191211_SG_4dFlow.dat'                                         'Sola_1',               'VE', 100; %31;
        'meas_MID00023_FID10176_BEAT_FQ_2148_4DFlowSG_PTEnabled_PT_on.dat'                          'PTPhan',               'VE', 100; %32;
        'meas_MID00069_FID10219_BEAT_FQ_2148_4DFlowSG_PTEnabled_PT_on.dat'                          'OPS',               'VE', 100; %33;
        'meas_MID00229_FID109571_BEAT_FQ_sg_whole_heart_8min.dat'                                   'NCH_FeraHeme_P2',     'VE', 100; %34;
        'meas_MID00057_FID24969_4D_flow_MV_PilotTone.dat'                                   'Sola_Mitral',     'VE', 100; %35;
        'meas_MID00487_FID23020_flow_RT_AO_PT_ON_notrig_match_sp_res_26s.dat'                                   'RT-PC-PT',     'VE', 100; %36;
        'meas_MID00099_FID34364_BEAT_FQ_2148_4DFlow_sg.dat'                                   'Sola_Pat1',     'VE', 100; %37;
        'meas_MID00301_FID43384_BEAT_FQ_2148_4DFlow_sg.dat'                                   'Sola_Pat2',     'VE', 100; %38;
        'meas_MID00399_FID143334_BEAT_FQ_sg_whole_heart_8min.dat'                                   'NCH_FeraHeme_P3',     'VE', 100; %39;
        'meas_MID00078_FID143013_BEAT_FQ_sg_whole_heart_8min.dat'                                   'NCH_FeraHeme_P4',     'VE', 100; %40;
        'meas_MID00081_FID144904_BEAT_FQ_sg_whole_heart_8min.dat'                                   'NCH_FeraHeme_P5',     'VE', 100; %41;
        'meas_MID00231_FID45443_BEAT_FQ_2148_4DFlow_sg.dat'                                   'Sola_Pat4',     'VE', 100; %42;
        'meas_MID00280_FID45492_BEAT_FQ_2148_4DFlow_sg.dat'                                   'Sola_Pat5',     'VE', 100; %43;
        'meas_MID00081_FID144904_BEAT_FQ_sg_whole_heart_8min.dat'                             'NCH_FeraHeme_P5',     'VE', 100; %44;
       
        
     
        

     };
 
saveName           = '4DSG';
% saveLocation       = 'D:\MRI Data\Aaron\4DSG MRM\binned\Cine RA Resp';
% saveLocation       = 'D:\MRI Data\Aaron\Joint Cine Flow (JCF)\1.5T Avanto\binned';
% saveLocation       = 'D:\MRI Data\Aaron\CardioStepper4D (CStep4D)\binned';
% saveLocation       = 'D:\MRI Data\Aaron\NCH FeraHeme 4D Flow (NCHFH4D)\binned';
saveLocation       = 'D:\MRI Data\Aaron\4DFlow SG vs PT\binned';
savePermutation    = 'Eff50';

fNames             = X(:,1);
tags               = X(:,2);
versions           = X(:,3);
resp_ranges        = X(:,4);

saveFigs = 0; 
% Dataset(s) to run
% *****************
ind = [9];
% *****************
%% =================================================================================================
% Options
% ==================================================================================================
opt = init();
% knobs
% ************************
opt.sigExt      = 'PCA';        % Signal extraction method, 'PCA' , 'SSAFARI', 'PT', or 'PTSola'
opt.ECG         = 0;
opt.rSoft       = 1;            % Soft-gating/weighting or not on respiratory dim
opt.nRPhases    = 1;            % # of respiratory phases to bin
opt.nCPhases    = 20;           % base # of cardiac phases to bin
opt.viewSharing = 0;            % view sharing to increase # of cardiac phases (x2)
opt.recStart    = 'start';     % 'start', 'center', 'end': where to begin trimming for "recTime"
opt.expOrd      = 4;
opt.respEff     = 50;
opt.flow        = 1;            % Flow dataset or not
opt.rampSS      = 25;           % Trim # self-gating samples during steady state ramp
% opt.rLPF        = [1];          % FeraHeme (P1)
% opt.cBPF        = [1.35 1.55];  % FeraHeme (P1)
% opt.rLPF        = [1.5];          % FeraHeme (Monkey)
% opt.cBPF        = [1.5 4.5];      % FeraHeme (Monkey)
% opt.rLPF        = [0.75];         % CardioStepper
% opt.cBPF        = [1.7 3];        % CardioStepper
% opt.cBPF        = [1.0 2.0]; 

% *************************
% ==================================================================================================
% Batch process 4D flow data
% ==================================================================================================
recTime = [300];
for i = 1 : length(ind)
     opt.fName = fNames{ind(i)};
     opt.tag   = tags{ind(i)};
     opt.version = versions{ind(i)};
     opt.resp_range = resp_ranges{ind(i)};
     
     cd(fileparts(which(opt.fName)));
     
     % Load raw data file and extract useful parameters
     disp('Loading raw data...')
     if ~strcmp(opt.sigExt, 'PTSola')
     [Dat, ~, param, ~, rawData] = readWrapperVBVDVE_SG2(opt, opt.version);
     else
     [Dat, ~, param, ~, rawData, pilot_tone_raw] = readWrapper_PTSola(opt);
     end
         
     
     
     
     for j = 1 : length(recTime)
        opt.recTime = recTime(j);
        
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % At some point need to create self-contained header
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        % Bin using self-gating signal only
        if ~strcmp(opt.sigExt, 'PTSola')
        FlowSG4D_Outputs_SG = binSelfGating(Dat, rawData, param, opt);
        else
        % Bin using self-gating
        opt.sigExt = 'PCA';        
        FlowSG4D_Outputs_SG = binSelfGating(Dat, rawData, param, opt);
        % Bin using pilot tone
        opt.sigExt = 'PTSola';
        FlowSG4D_Outputs_PT = binSelfGating_PTSola(Dat, rawData, param, opt, pilot_tone_raw);
        end
%         Signals.SG.Time = FlowSG4D_Outputs_SG.Time + opt.recTime/2;
%         Signals.SG.Respiratory = FlowSG4D_Outputs_SG.Respiratory;
%         Signals.SG.Cardiac = FlowSG4D_Outputs_SG.Cardiac;
%         Signals.PT.Time = FlowSG4D_Outputs_PT.Time + opt.recTime/2;
%         Signals.PT.Respiratory = FlowSG4D_Outputs_PT.Respiratory;
%         Signals.PT.Cardiac = FlowSG4D_Outputs_PT.Cardiac;       
        
        
        % Construct array (2)
        name = [saveName,'_',savePermutation,'_',opt.tag,'_',num2str(opt.recTime/60),'min'];
%         name = [saveName,'_',savePermutation,'_',opt.tag,'_',num2str(opt.recTime),'sec'];
        raw_header = rawData.hdr;
        [resCol,resLin,resPar,rate_samps,rate_weights,NImageCols,NImageLins,NImagePars] = constructArrayWrapper(FlowSG4D_Outputs_SG, opt, name, saveLocation, param, raw_header);
        
%         FlowSG4D_Outputs.param.ResCol = resCol;
%         FlowSG4D_Outputs.param.ResLin = resLin;
%         FlowSG4D_Outputs.param.ResPar = resPar;
%         FlowSG4D_Outputs.param.rate_samps = rate_samps;
%         FlowSG4D_Outputs.param.rate_weights = rate_weights;
%         FlowSG4D_Outputs.param.NImageCols = NImageCols;
%         FlowSG4D_Outputs.param.NImageLins = NImageLins;
%         FlowSG4D_Outputs.param.NImagePars = NImagePars;       

        % Save results
%         name = [saveName,'_',opt.tag,'_',num2str(opt.recTime/60),'min'];
% %         name = [saveName,'_',opt.tag,'_',num2str(opt.recTime),'sec'];
%         save(name,'FlowSG4D_Outputs','-v7.3');
        
        % Save figures        
        if saveFigs
            h =  findobj('type','figure');
            n = length(h);
            
            if ~exist([saveLocation,'\',name,'\figs'],'dir')
                mkdir([saveLocation,'\',name,'\figs']);
            end
            pause();
            for fig = 1 : n
                filename = [saveLocation,'\',name,'\figs\','Fig_',num2str(fig)];
                frame = get(handle(h(fig)),'JavaFrame');
                exportFigure(filename, frame);
            end
        end
     end
     
end





function exportFigure(fileName, frame_h)
pause(0.00001);
% frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
snam='default'; % The name of your style file (NO extension)
s=hgexport('readstyle',snam);
%apply style sheet info
fnam=fileName; % your file name
s.Format = 'tiff'; %I needed this to make it work but maybe you wont.
hgexport(gcf,fnam,s);
pause(0.00001)
close(gcf);
end