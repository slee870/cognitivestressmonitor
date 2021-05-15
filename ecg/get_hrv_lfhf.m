function [lfhf,tref,lfr,hfr] = get_hrv_lfhf(ecg,fs,twin,overlap)
%%% * Compute heart rate variability based on low-freq to high-freq ratio,
%%% as done in Vanoli 1995
% IN
%     - ecg : Nx1 vector of raw ECG
%     - fs : sampling rate (hz)
%     - twin : window size for calculation (s)
%     - overlap : overlap between sliding windows ([0 -> 1])
% OUT
%     - lfhf : Mx1 vector of low-freq / high-freq ratio
%     - tref : time index of frames within which lfhf is calculated
%     - lfr : low-freq ratio (% of full band)
%     - hfr : high-freq ratio (% of full band)

detrendOn = 0;
lf_band     = [0.04 0.15];
hf_band     = [0.15 5];
full_band   = [0.04 5];

% EXTRACT AC AND DC COMPONENTS
if detrendOn == 1
    % dcwin = 10;
    dcwin   = twin*5;
    ecg     = extract_acdc(ecg,dcwin,fs);
end

% DIVIDE SIGNALS INTO INTERVAL FRAMES
len             = length(ecg);
framesize       = round(twin * fs);
shift           = round(framesize * (1-overlap));

i1      = 1;
incr    = 1;
while i1 < (len - framesize)
    i2              = i1 + framesize;
    frame           = i1:i2;
    ecg_snippet     = ecg(frame);
    
    lf = bandpower(ecg_snippet,fs,lf_band);
    hf = bandpower(ecg_snippet,fs,hf_band);
    full = bandpower(ecg_snippet,fs,full_band);
    
    lfhf(incr,1) = lf ./ hf;
    lfr(incr,1) = lf ./ full;
    hfr(incr,1) = hf ./ full;
    
    iref(incr)  = i1;
    tref(incr)  = i1 / fs;
    i1          = i1 + shift;
    incr        = incr + 1;
end

% while i1 >= (len - framesize) && i1 < len
%     
% end

tref  = tref';

end