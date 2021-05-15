function [iref,tref] = detect_R_peaks_threshold(ecg_f,fs,pct,plotOn)
% *** Detect R-wave peaks on ECG based on fixed threshold envelope
% IN
%     - ecg_f : Nx1 vector of filtered ECG signal
%     - fs : sampling rate of ecg_f (Hz)
%     - pct : threshold for peak detection as pct of signal peak envelope
%         0 -> 1
%     - plotOn : plot result (1) or not (0)
% OUT
%     - iref : vector of R peak indices (samples)
%     - tref : vector of R peak times (s)
    
win = 5*fs;
dur_min = 0.4;  % s
if length(ecg_f) < win
    win = length(ecg_f);
end

env = envelope_nbb(ecg_f,win,'peak');
minpeakdist = round(dur_min*fs);
minpeakheight = pct * median(env);
% minpeakprom = ;

% DETECT PEAKS ON LOGICAL SIGNAL
% on_logic = ecg_f > minpeakheight;
% on = double(on_logic);
% [px,ix] = findpeaks(on,'minpeakdistance',minpeakdist);

% DETECT PEAKS ON ECG SIGNAL ITSELF
[~,ix] = findpeaks(ecg_f,'minpeakdistance',minpeakdist,'minpeakheight',minpeakheight);

iref = ix;
tref = iref / fs;

% if plotOn
%     figure;
%     sp(1) = subplot(311);
%     hold on;
%     plot((1:length(ecg_f))/fs,ecg_f)
%     plot(tref,ecg_f(ix),'.','markersize',10);
%     legend('ecg','R peak');
%     xlabel('Time (s)');
%     sp(2) = subplot(312);
%     hold on;
% %     plot(tref(1:end-1),diff(tx),'.','markersize',10);
%     plot(tref(1:end-1),diff(tref),'.','markersize',10);
%     sp(3) = subplot(313);
%     plot(tref(1:end-2),diff(diff(tref)),'.','markersize',10);
%     linkaxes(sp,'x');
% end

if plotOn
    figure;
    
   % sp(1) = subplot(211);
    hold on;
    plot((1:length(ecg_f))/fs,ecg_f, 'Linewidth' ,2)
    plot(tref,ecg_f(iref),'.','markersize',30);
    legend('ecg','R peak');
    xlabel('Time (s)');
    ylabel('Normalized AU')
    title('R Peak Detection');
    box on
    hhh = 60 ./ diff(tref);
% h = smoothdata(hhh,'Gaussian', 5);
% sp(2) = subplot(212);
% plot(tref(1:end-1),h+20,'b','Linewidth', 2 );
  %  xlabel('Time (s)');
  %  ylabel('HR (BPM)')
  %  title('HR');
    
 %   linkaxes(sp,'x');
    set(gcf,'color','w')
    box on
end
end