close all
warning off
clear

dat_fold = 'csvs\apnea\'; % Data folder
folders = {'kl_3_25_21', 'kl_4_8_21', 'nz_3_25_21', 'nz_4_8_21', 'shl_3_25_21', 'shl_4_8_21' }; % subject folders
flabels = {'mecg_', 'ppg_'};
ppg_fs = 50.0; acc_fs = 250.0; % Sampling rates
FsNew = 5;

%% Loop
for s = 1%:length(folders)
    %% Load Data
    T1 = readtable([flabels{1}, folders{s} '.csv']); % load motion/ecg data
    T2 = readtable([flabels{2}, folders{s} '.csv']); % load ppg data
    
    accdata = table2array(T1(:,1:5)); %Break into arrays
    ppgdata  = table2array(T2);
    
    tp =  ppgdata(:,1); %Time for PPG
    ch1 = ppgdata(:,2); %IR
    ch2 = ppgdata(:,3); %Red
    tm =  accdata(:,1); %Time for motion-ecg
    ax =  accdata(:,2);
    ay =  accdata(:,3);
    az =  accdata(:,4);
    ecg = accdata(:,5);
    
    e = [tm, ecg]; %Concatendate ECG and time
    ind = ecg < 650; %Remove placeholders
    e = e(ind,:); %Remove placeholders
    
    %% Upsample
    ind = find(isfinite(ch1));
    ch1= ch1(ind);
    ch2 = ch2(ind);
    ind = isnan(ch1);
    ch1 = ch1(~ind);
    ch2 = ch2(~ind);
    tp = tp(~ind);
    ind = isnan(ch2);
    ch1 = ch1(~ind);
    ch2 = ch2(~ind);
    tp = tp(~ind);
    
    ind = find(isfinite(az));
    az = az(ind);
    ay = ay(ind);
    ax = ax(ind);
    tm = tm(ind);
    ind = isnan(az);
    az = az(~ind);
    az = az(~ind);
    ax = ax(~ind);
    tm = tm(~ind);
    
    ir = upsample_ppg(ch1, tp, tm); %resample everything at 250
    red = upsample_ppg(ch2, tp, tm);
    ecg2 = upsample_ppg(e(:,2), e(:,1), tm);
    
    ind = find(isfinite(ir));
    ir= ir(ind);
    red = red(ind);
    ind = isnan(ir);
    ir = ir(~ind);
    red = red(~ind);
    
    
    %% Filter
    
    cutoff = [.8 5]; %BPF for PPG
    Fs = 250;
    [bh, ah] = butter(1, cutoff*2/Fs, 'bandpass');
    red_f = filtfilt(bh, ah, red);
    ir_f = filtfilt(bh, ah, ir);
    
    cutoff = [4 20]; %BPF for SCG
    [bh, ah] = butter(1, cutoff*2/Fs, 'bandpass');
    ac_f = filtfilt(bh, ah, az);
    
    %resperation
    cutoff = [0.01 0.4]; %BPF for resp
    Fs = 250;
    [bh, ah] = butter(1, cutoff*2/Fs, 'bandpass');
    resp = filtfilt(bh, ah, az);
    
    %% Spo2
    
    [spo2s, spo2t] = rmsSpO2_v2(ch1, ch2, 50);
    
    %% ECG
    
    [HR, HRV, ecgn, thrv, rix, ecgf] = ECG_params(ecg2, tm, 250); %Get HR and HRV
    [rix]   = detect_R_peaks_threshold(ecgn,Fs,0.7,1); %Find each R peak from ECG
    tref    = rix / Fs; %R peaks in actual time
    tHR = 1:length(HR);
    
           %% SCG
    
        [~,scg_mat] = segment_cycles(ac_f,tref,Fs,0); %Segment each SCG by the R peak
    
        plotOn = 1;
        twin = 6;
        overlap = 1/3;
        [scg_avg,tref_avg,beats_used] = avg_beat_interval(ac_f,tref,Fs,twin,overlap,plotOn); %segmented SCG based on beat averages
    
        SCG = scg_avg;
        s1win = [0.1 0.2];
        s2win = [0.45 0.6];
        lbl = label_scg(SCG,Fs,s1win,s2win,1); %picking out fiducials - this is where the magic happens
        PEP = lbl.ao;
        LVET = lbl.ir - lbl.ao;
        PEP_LVET = PEP ./ LVET;
        AO_magn = lbl.aom;
        S1_ht = lbl.aom - lbl.icm;
        S1_S2 = lbl.aom ./ lbl.ir;
        
%          
% av1 = mean(SCG',2,'omitnan');
% 
%  
% sd1 = std(SCG,'omitnan');
% t = (1:299)./50;
% 
% figure
% plot(t,av1,'r');
% hold on
% plotFillAlpha(t,av1,sd1,'r',0.2);
    
    %% Respiration
    [RR, RR_t, filt] = findRR (ac_f, tm, 250);
    plot(RR_t, RR)
    
    % Signal quality metrics
    [cPPG,tref_cPPG] = get_moving_cardiodicity(red_f,250,4,0.5);
    [cResp,tref_cResp] = get_moving_cardiodicity(resp,250,4,0.5);
    [cSCG,tref_cSCG] = get_moving_cardiodicity(ac_f,250,4,0.5);
    [yup1,ylow1] = envelope(ac_f,500, 'rms');
    figure; plot(yup1);
    [yup2,ylow2] = envelope(az,500, 'rms');
     figure; plot(yup2);
    [yup3,ylow3] = envelope(az,500, 'peak');
    figure; plot(yup3);
    [yup4,ylow4] = envelope(az,500, 'analytic'); %analytic reconstruction in orthoganol hilbert space
    figure; plot(yup4);
    title('Hilbert Analytic Envelope of Acceleration')
    xlabel('Time (s)')
    ylabel('Acceleration [uG]')
    set(gcf,'color','w');
    [yup5,ylow4] = envelope(az,500, 'analytic'); %analytic reconstruction in orthoganol hilbert space
    
    % resmaple everything at 5Hz
    
    new_t = 1:0.2:max(tm);
    yup1_rs = interp1(tm, yup1, new_t, 'spline');
    yup2_rs = interp1(tm, yup2, new_t, 'spline');
    yup3_rs = interp1(tm, yup3, new_t, 'spline');
    yup4_rs = interp1(tm, yup4, new_t, 'spline');
    yup5_rs = interp1(tm, yup5, new_t, 'spline');
    cSCG_rs = interp1(tref_cSCG, cSCG, new_t, 'spline');
    cPPG_rs = interp1(tref_cPPG, cSCG, new_t, 'spline');
    RR_rs = interp1(RR_t, RR, new_t, 'spline');
    HRV_rs = interp1(thrv, HRV, new_t, 'spline');
    HR_rs = interp1(tHR, HR, new_t, 'spline');
    
    %% save all this junk
    
    fname = [folders{s}, '_features.mat'];
    save(fname, 'yup1', 'yup2', 'yup3', 'yup4', 'yup5','cSCG_rs', 'cPPG_rs', 'RR_rs', 'HRV_rs', 'HR_rs')
end






