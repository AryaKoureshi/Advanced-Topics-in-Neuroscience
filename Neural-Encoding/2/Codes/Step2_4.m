clc;clear;close all;

%% Loading data and initializing
load 'UnitsData.mat';
NumBins = 64; % number of bins
WidthWindow = 9; % width of the window
cntrs = linspace(-1.2, 2, NumBins);
ShuffleStatus = input('Set status of the shuffling: (1=on, 0=off) '); % shuffle the final labels
NeuronIdxs = 1:numel(Unit);
pValues = [];
VecIdx = 1:numel(Unit)*2;
VecIdx = VecIdx(randperm(length(VecIdx)));
y = [];
CntAll = [];

%% Step 2_4 - Population - Modeling LR cue
for NeuronIdx = NeuronIdxs
    IdxTrials1 = [];
    IdxTrials2 = [];
    
    pos = -1;
    for i = 1:numel(Unit(NeuronIdx).Cnd)
        cnd = Unit(NeuronIdx).Cnd(i);
        if cnd.Value(2) == pos
            IdxTrials1 = [IdxTrials1; cnd.TrialIdx];
        end
    end
    data = Unit(NeuronIdx).Trls(IdxTrials1);
    [cnts,~] = PSTH(data, WidthWindow, NumBins, cntrs);
    CntAll = [CntAll; cnts];
    y = [y; 0];

    pos = 1;
    for i = 1:numel(Unit(NeuronIdx).Cnd)
        cnd = Unit(NeuronIdx).Cnd(i);
        if cnd.Value(2) == pos
            IdxTrials2 = [IdxTrials2; cnd.TrialIdx];
        end
    end
    data = Unit(NeuronIdx).Trls(IdxTrials2);
    [cnts,~] = PSTH(data, WidthWindow, NumBins, cntrs);
    CntAll = [CntAll; cnts];
    y = [y; 1];
end

if ShuffleStatus
    y = y(VecIdx);
end

mdl = fitglm(CntAll, y);
LR_pValuePopulation = coefTest(mdl)

%% Functions
function [cnts, cntrs] = PSTH(data, WidthWindow, NumBins, cntrs)
    data_all = zeros(numel(data), NumBins);
    for i=1:numel(data)
        [cnts, cntrs] = hist(cell2mat(data(i)), NumBins, 'xbins', cntrs);
        cnts = movmean(cnts, WidthWindow);
        data_all(i, :) = cnts;
    end
    data_all = mean(data_all,1);
    cnts = data_all/(3.2/NumBins);

end

function plot_PSTH(Unit, NeuronIdx, WidthWindow, NumBins, cntrs)
    hold on
    CntAll = [];
    for CueValue = 3:3:9
        for pos = [-1, 1]
            value = [CueValue, pos];
            data = GetCnd(Unit, NeuronIdx, value);
            [cnts,~] = PSTH(data, WidthWindow, NumBins, cntrs);
            CntAll = [CntAll; cnts];
            plot(cntrs, cnts, 'LineWidth', 1)
        end
    end
    CntAll = mean(CntAll,1);
    plot(cntrs, CntAll, 'k', 'LineWidth', 1.5)
    title("PSTH of Unit " + num2str(NeuronIdx))
    xlabel("Time (sec)")
    ylabel('Firing Rate (Hz)')
    xlim([-1.2, 2])
    xline(0, '--', 'Reward Cue', LabelHorizontalAlignment = 'center');
    xline(0.3, '--', 'Delay Period', LabelHorizontalAlignment = 'center');
    xline(0.9, '--', 'Reaction', LabelHorizontalAlignment = 'center');
    legend('[3 -1]', '[3 1]', '[6 -1]', '[6 1]', '[9 -1]', '[9 1]', 'Average')
    hold off
end

function data = GetCnd(Unit, NeuronIdx, value)
    for i = 1:numel(Unit(NeuronIdx).Cnd)
        cnd = Unit(NeuronIdx).Cnd(i);
        if cnd.Value == value
            trials_indx = cnd.TrialIdx;
            data = Unit(NeuronIdx).Trls(trials_indx);
            break
        end
    end
end