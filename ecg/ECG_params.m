function [HR, HRV,ecg_filtered_norm, thrv, Rpks, ecg_filtered ] = ECG_params(ecgData, time, Fs)

% Need to first filter the signals to isolate QRS complex
cutoff = [5 40];
[bh, ah] = butter(1, cutoff*2/Fs, 'bandpass');
ecg_filtered = filtfilt(bh, ah, ecgData);

% ecg_filtered = ecg_filtered.*(-1); %If ECG leads are backwards :0

% Normalize the ECG
ecg_filtered_norm=ecg_filtered-mean(ecg_filtered);
ecg_filtered_norm=ecg_filtered_norm./movmax(ecg_filtered_norm, Fs*3);
ecg_filtered_norm(ecg_filtered_norm > 10) = 0;
ecg_filtered_norm(ecg_filtered_norm < -10) = 0;
ecg_filtered_norm(isnan(ecg_filtered_norm)) = 0;

min_peak_height=0.5;

t = 1:max(time);
% Finds the peaks
[pks,ecg_pks]=findpeaks(ecg_filtered_norm,'MinPeakHeight',min_peak_height,'MinPeakDistance',150);

Rpks = time(ecg_pks);

ibi = diff(Rpks);
k = 3;
if k == 1
    HR = 60 ./ ibi;
else
    HR = 60 ./ movmean(ibi,k);
end

% figure
% plot(HR)

[HRV,thrv] = get_hrv(Rpks,5,0.5,'rmssd',1);

end