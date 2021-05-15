function [iref,tref] = detect_R_peaks_threshold_sweep(ecg_f,fs,method,plotOn)
% *** Detect R-wave peaks on ECG based on a swept threshold envelope (i.e.,
% choose the optimal threshold based on 'method' criterion)
% IN
%     - ecg_f : Nx1 vector of filtered ECG signal
%     - fs : sampling rate of ecg_f (Hz)
%     - method : by what metric should R peak detection be judged?
%         'corrcoef' -> max R^2 (i.e., max linearity of tref)
%         'std' -> min st. dev. (i.e., least variation beat-to-beat in tref)
%     - plotOn : plot result (1) or not (0)
% OUT
%     - iref : vector of R peak indices (samples)
%     - tref : vector of R peak times (s)

win = 5*fs;
dur_min = 0.3;  % s
if length(ecg_f) < win
    win = length(ecg_f);
end

env = envelope_nbb(ecg_f,win,'peak');
minpeakdist = round(dur_min*fs);
% minpeakprom = ;

pcts = [0.5:0.05:0.95]; % sweep envelope threshold from 50% to 95%
i_pct = 1:length(pcts);
for i = i_pct
    minpeakheight = pcts(i) * median(env);
    
    % DETECT PEAKS ON ECG SIGNAL ITSELF
    [~,ix{i}] = findpeaks(ecg_f,'minpeakdistance',minpeakdist,'minpeakheight',minpeakheight);
    
    % DETECT PEAKS ON LOGICAL SIGNAL
    % on_logic = ecg_f > minpeakheight;
    % on = double(on_logic);
    % [px,ix] = findpeaks(on,'minpeakdistance',minpeakdist);

    switch method
        case 'corrcoef'
            coeffs = corrcoef(1:length(ix{i}),ix{i});
            crit(i) = coeffs(2);
        case 'std'
            crit(i) = std(diff(ix{i}));
    end
end

switch method
    case 'corrcoef'
        [opt,iopt] = max(crit);
    case 'std'
        [opt,iopt] = min(crit);
end

iref = ix{iopt};
tref = iref / fs;

if plotOn
    figure;
    
    sp(1) = subplot(311);
    hold on;
    plot((1:length(ecg_f))/fs,ecg_f)
    plot(tref,ecg_f(iref),'.','markersize',10);
    legend('ecg','R peak');
    xlabel('Time (s)');
    title('Threshold = ' + string(pcts(iopt)));
    
    sp(2) = subplot(312);
    plot(tref(1:end-1),60 ./ diff(tref),'.','markersize',10);
    title('HR');
    
    subplot(313);
    histogram(60 ./ diff(tref));
    xlabel('Instantaneous HR');
    ylabel('Count');
    title(string(method) + ' = ' + string(opt));

    linkaxes(sp,'x');
end
end