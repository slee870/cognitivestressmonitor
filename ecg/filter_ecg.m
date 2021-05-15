function [ecg_f] = filter_ecg(ecg,fs,plotOn)

Fstop = 60;

% Lowpass IIR
iir_order = 7;
iir = designfilt('lowpassiir','filterorder',iir_order,...
    'halfpowerfrequency',Fstop,'samplerate',fs);
ecg_f = filtfilt(iir,ecg);

if plotOn
    figure;
    hold on;
    plot(ecg);
    plot(ecg_f);
end

end