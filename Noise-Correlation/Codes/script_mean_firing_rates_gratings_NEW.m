%  This script contains three parts:
%   1. Convert spike times to 1ms bins.
%   2. Compute the trial-averaged population activity (PSTHs).
%
%  S struct is used to store the data:
%    S(igrat).trial(itrial).spikes  (num_units x num_1ms_timebins)
%    S(igrat).trial(itrial).counts  (num_units x num_20ms_timebins)
%    S(igrat).mean_FRs  (num_units x num_20ms_timebins)
%
%  Author: Ben Cowley, bcowley@cs.cmu.edu, Oct. 2016
%
%  Notes:
%    - automatically saves 'S' in ./pvc-11/data_and_scripts/spikes_gratings/


%% parameters

    SNR_threshold = 1.5;
    firing_rate_threshold = 1.0;  % 1.0 spikes/sec
    binWidth = 20;  % 20 ms bin width

    
%% parameters relevant to experiment

    length_of_gratings = 1;  % each gratings was shown for 1.28s, take the last 1s
    
    filenames{1} = './pvc-11/data_and_scripts/spikes_gratings/data_monkey1_gratings.mat';
    filenames{2} = './pvc-11/data_and_scripts/spikes_gratings/data_monkey2_gratings.mat';
    filenames{3} = './pvc-11/data_and_scripts/spikes_gratings/data_monkey3_gratings.mat';


    monkeys = {'monkey1', 'monkey2', 'monkey3'};
    

%%  spike times --> 1ms bins

    for imonkey = 1:length(monkeys)
        S = [];

        fprintf('binning spikes for %s\n', monkeys{imonkey});

        load(filenames{imonkey});
            % returns data.EVENTS

        num_neurons = size(data.EVENTS,1);
        num_gratings = size(data.EVENTS,2);
        num_trials = size(data.EVENTS,3);

        edges = 0.28:0.001:1.28;  % take 1ms bins from 0.28s to 1.28s

        for igrat = 1:num_gratings
            for itrial = 1:num_trials
                for ineuron = 1:num_neurons
                    S(igrat).trial(itrial).spikes(ineuron,:) = histc(data.EVENTS{ineuron, igrat, itrial}, edges);
                end
                S(igrat).trial(itrial).spikes = S(igrat).trial(itrial).spikes(:,1:end-1);  % remove extraneous bin at the end
            end
        end

        save(sprintf('./pvc-11/data_and_scripts/spikes_gratings/S_%s.mat', monkeys{imonkey}), 'S', '-v7.3');
    end
   
%%  Take spike counts in bins
    for imonkey = 1:length(monkeys)
        
        fprintf('spike counts in %dms bins for %s\n', binWidth, monkeys{imonkey});
        
        load(sprintf('./pvc-11/data_and_scripts/spikes_gratings/S_%s.mat', monkeys{imonkey}));
            % returns S(igrat).trial(itrial).spikes
        num_grats = length(S);
        num_trials = length(S(1).trial);
        
        for igrat = 1:num_grats
            for itrial = 1:num_trials
                S(igrat).trial(itrial).counts = bin_spikes(S(igrat).trial(itrial).spikes, binWidth);
            end
        end
        
        save(sprintf('./pvc-11/data_and_scripts/spikes_gratings/S_%s.mat', monkeys{imonkey}), 'S', '-v7.3');
    end
    
    
%%  Compute trial-averaged population activity (PSTHs)

    for imonkey = 1:length(monkeys)
        fprintf('computing PSTHs for %s\n', monkeys{imonkey});
        
        load(sprintf('./pvc-11/data_and_scripts/spikes_gratings/S_%s.mat', monkeys{imonkey}));
            % returns S(igrat).trial(itrial).spikes
        num_grats = length(S);
        num_trials = length(S(1).trial);
        
        for igrat = 1:num_grats
            mean_FRs = zeros(size(S(igrat).trial(1).counts));
            for itrial = 1:num_trials
                mean_FRs = mean_FRs + S(igrat).trial(itrial).counts;
            end
            S(igrat).mean_FRs = mean_FRs / num_trials;
        end
        
        save(sprintf('./pvc-11/data_and_scripts/spikes_gratings/S_%s.mat', monkeys{imonkey}), 'S', '-v7.3');
    end
    

        

