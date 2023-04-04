function [S] = computeSpontCounts(EVENTS,bininseconds,spacebetweenbins)
% 
% [S, channels, medSNR] = computeSpontCounts(EVENTS,bininseconds,chunksize)
%
% Computes spike counts from spiking event times. 
%
% Inputs: 
%     EVENTS:           Cell array of spiking events for each sorted unit
%     bininseconds:     Number of seconds over which spike counts will be
%                       binned
%     spacebetweenbins: Time between successive bins.
%
% Outputs:
%     S:                N X T matrix of spike counts with N neurons and T 
%                       trials. Adjacent columns represent spike counts for
%                       adjacent time bins of size bininseconds
% Author: Ryan Williamson, Oct. 2016
%

 
maxtime = max(cell2mat(cellfun(@(x)max(x),EVENTS,'UniformOutput',false)));
maxtrials = floor(maxtime/(bininseconds+spacebetweenbins));
S =zeros(length(EVENTS), maxtrials);
for n = 1:maxtrials
    starttrial = (n-1)*(bininseconds+spacebetweenbins);
    endtrial = starttrial+bininseconds;
    S(:,n) = cell2mat(cellfun(@(x)sum(x>=starttrial&x<endtrial),EVENTS,'UniformOutput',false))';
end

end
