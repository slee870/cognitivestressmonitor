function [bpm] = get_hr(Rpks,k)
% *** Calculate heart rate (in bpm) from inter-beat interval
% IN
%     - Rpks : 1xN vector of R-peak timing (s) 
%         N -> # of beats
%     - k : # beats over which to average
% OUT
%     - bpm : heart rate (bpm)

if nargin < 2
    k = 1;
end

ibi = diff(Rpks);
if k == 1
    bpm = 60 ./ ibi;
else
    bpm = 60 ./ movmean(ibi,k);
end

bpm(end+1) = bpm(end);  % just repeat the last to shore up numbers

end