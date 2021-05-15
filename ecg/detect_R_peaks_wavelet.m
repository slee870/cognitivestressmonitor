function [Rpx,Rix] = detect_R_peaks_wavelet(ecg,fs,plotOn)
% *** Extract R peaks from ECG data for beat segmentation
% IN
%     - tbl : table containing at least time and ECG (tbl.time, tbl.ecg)
%     - fs : sampling rate (Hz)
%     - plotOn : (1) plot output or (0) don't
% OUT
%     - Rpks : R-wave peaks [s]
%     - grpdel : group delay introduced by filter

%% FILTER ECG
% FILTER METHOD 1
ecg_f = filter_ecg(ecg,fs,0);

% FILTER METHOD 2
% ecg_band = [0 0.5 50 60];
% [ecg_f,grpdel] = apply_filter(ecg,ecg_band,fs,'bpf');
% ecg_shift = ecg_f(1+grpdel_end));

%% SIMPLE PEAK DETECTION
% minpeakprom = 2*rms(ecg_f);
% minpeakdist = 0.25*fs;
% [Rpx,Rix] = findpeaks(ecg_f,'minpeakprominence',minpeakprom,'minpeakdistance',minpeakdist);

%% WAVELET-TRANSFORMED ECG TO ISOLATE QRS COMPLEX
ecg_w = modwt(ecg_f,5);             % decomp. ECG down to level 5
ecg_w_rc = zeros(size(ecg_w));      % reconstruct ECG at scales 4 and 5
ecg_w_rc(4:5,:) = ecg_w(4:5,:);
qrs = imodwt(ecg_w_rc,'sym4');
qrs_2 = abs(qrs).^2;

% DETECT PEAKS FROM ISOLATED QRS
minpeakprom = quantile(abs(qrs),0.995);
minpeakdist = 0.35*fs;
[Rpx,Rix] = findpeaks(qrs,'minpeakprominence',minpeakprom,'minpeakdistance',minpeakdist);

if plotOn
    figure, hold on
    plot(ecg)
    plot(ecg_f)
    plot(qrs)
    plot(Rix,ecg_f(Rix),'o')
    legend('raw','filt','qrs');
end


