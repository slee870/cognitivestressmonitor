function [hrv,tref] = get_hrv(Rpks,win,ovr,method,plotOn)
%%% * Compute heart rate variability given R peaks (or other reliable
%%% cardiac timing interval)
% IN
%     - Rpks : Nx1 vector of cardiac timing indices (s)
%     - win : window over which to average or otherwise compute HRV (s)
%         * win == [] -> defaults to full range of input
%     - ovr : overlap between successive windows (pct, [0 1])
%         * ovr == [] -> defaults to 0 %
%     - method : method by which to compute HRV
%         * 'rmssd' -> root-mean-square of successive differences 
%         * 'sdnn' -> standard deviation of beat intervals 
%     - plotOn : plot result (1) or not (0)
% OUT
%     - hrv : heart rate variability metric within given frame
%     - tref : time index of frames within which HRV is calculated
    
if isrow(Rpks)
    Rpks = Rpks';
end 

if isempty(win)
    win = range(Rpks);
end

if isempty(ovr)
    ovr = 0;
end
   
switch method
    case 'rmssd'    % root mean square of successive IBI differences
        ibi     = diff(Rpks);   % inter-beat interval
        ibi     = [ibi; ibi(end)];
        sd      = diff(ibi);     % successive differences
        sd      = [sd; sd(end)];
        ssd     = sd.^2;        % square of successive differences        
        [rmssd,tref] = win_rmssd(Rpks,ssd,win,ovr); % compute RMSSD
        hrv     = rmssd;
    case 'sdnn'     % standard deviation of normal-normal IBI
        ibi     = diff(Rpks);   % inter-beat interval
        ibi     = [ibi; ibi(end)];
        [sdnn,tref] = win_sdnn(Rpks,ibi,win,ovr);   % compute SDNN
        hrv = sdnn;
end

if ~exist('tref','var')
    tref = [];
end

if plotOn
    figure;
    hold on;
    sp(1) = subplot(211);
    plot(Rpks, 60 ./ ibi);
    title('Heart Rate (bpm)');
    sp(2) = subplot(212);    
    plot(tref,hrv);
    title(['HRV (' method ')']);
    linkaxes(sp,'x');
end 
    
    
end

%% HELPERS
%%% RMSSD WINDOWING
function [rmssd,tref] = win_rmssd(Rpks,ssd,win,ovr)

t1 = Rpks(1);
tshift = win * (1-ovr);
incr = 1;

while t1 < (Rpks(end) - win)    
    t2 = t1 + win;
    ix_inframe = Rpks > t1 & Rpks < t2;
    ssd_inframe = ssd(ix_inframe);    
    mssd = mean(ssd_inframe);
    rmssd(incr,1) = sqrt(mssd);
    tref(incr,1) = t1;    
    t1 = t1 + tshift;
    incr = incr + 1;
end

% END CONDITION
if t1 > (Rpks(end) - win)
    ix_inframe = Rpks > t1 & Rpks <= Rpks(end);
    ssd_inframe = ssd(ix_inframe);    
    mssd = mean(ssd_inframe);
    rmssd(incr,1) = sqrt(mssd);
    tref(incr,1) = t1;
end 

end

%%% SDNN WINDOWING
function [sdnn,tref] = win_sdnn(Rpks,ibi,win,ovr)

t1 = Rpks(1);
tshift = win * (1-ovr);
incr = 1;

while t1 < (Rpks(end) - win)    
    t2 = t1 + win;
    ix_inframe = Rpks > t1 & Rpks < t2;
    ibi_inframe = ibi(ix_inframe);    
    sdnn(incr,1) = std(ibi_inframe);
    tref(incr,1) = t1;    
    t1 = t1 + tshift;
    incr = incr + 1;
end

% END CONDITION
if t1 > (Rpks(end) - win)
    ix_inframe = Rpks > t1 & Rpks <= Rpks(end);
    ibi_inframe = ibi(ix_inframe);    
    sdnn(incr,1) = std(ibi_inframe);
    tref(incr,1) = t1;
end 

end
