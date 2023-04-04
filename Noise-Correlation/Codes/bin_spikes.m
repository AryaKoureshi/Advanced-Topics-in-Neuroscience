function spikecounts = bin_spikes(spikes, binWidth)
% S = bin_spikes(spikes, binWdith)
% 
% bins the spikes with given binWidth (in ms)
%
% INPUT:
%   spikes (num_neurons x num_1ms_timepoints):  contains binary
%       spike information in 1ms bins
%   binWidth (in ms): bins spikes in non-overlapping windows with
%       specified bin width
%
% OUTPUT:
%   spikecounts (num_neurons x num_timebins): contains
%       binned spike counts
%
% NOTES:
%   If a trial has a number of timepoints not evenly divisible by
%       the binWidth, the remainder timepoints are removed and not
%       considered during the spiking.
%
%
% Author: bcowley, 2014


    spikecounts = [];
    
    num_timepoints = size(spikes,2);
    
    if (num_timepoints < binWidth)
        warning('bin_spikes:  Number of timepoints is less than binWidth.');
    end
            
    bin_indices = 1:binWidth:num_timepoints;

    if (mod(num_timepoints,binWidth) ~= 0 && length(bin_indices) > 1) % the last bin does not have enough time points in it
        bin_indices = bin_indices(1:(end-1));
    end

    for ibin = 1:length(bin_indices)
        
        spikecounts = [spikecounts ...
            sum(spikes(:, bin_indices(ibin):(bin_indices(ibin)+binWidth -1)),2)];
    end            
 

end



