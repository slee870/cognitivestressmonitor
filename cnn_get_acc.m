close all
warning off
clear

dat_fold = 'csvs\apnea\'; % Data folder
folders = {'kl_3_25_21', 'kl_4_8_21', 'nz_3_25_21', 'nz_4_8_21', 'shl_3_25_21', 'shl_4_8_21' }; % subject folders
flabels = {'mecg_', 'ppg_'};
ppg_fs = 50.0; acc_fs = 250.0; % Sampling rates
FsNew = 5;

%% Loop
for s = 1:length(folders)
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
    
    cutoff = [0.1 50]; %BPF for SCG
    Fs = 250;
    [bh, ah] = butter(1, cutoff*2/Fs, 'bandpass');
    ac_f = filtfilt(bh, ah, az);
        
    fname = [folders{s}, '_acc.mat'];
    save(fname, 'ac_f')
end
