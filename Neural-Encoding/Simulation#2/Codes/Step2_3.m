clc;clear;close all;

%% Loading data and initializing
load 'UnitsData.mat';
NumBins = 64; % number of bins
WidthWindow = 9; % width of the window
cntrs = linspace(-1.2, 2, NumBins);
ShuffleStatus = input('Set status of the shuffling: (1=on, 0=off) '); % shuffle the final labels
NeuronIdxs = 1:numel(Unit);
pValues = [];
VecIdx = 1:numel(Unit(1).Trls);
VecIdx = VecIdx(randperm(length(VecIdx)));

%% Step 2_3 - Single Units - Modelling LR cue and reward expected value
for NeuronIdx = NeuronIdxs
    CntAll = [];
    IdxTrials1 = [];
    IdxTrials2 = [];
    IdxTrials3 = [];
    IdxTrials4 = [];
    IdxTrials5 = [];
    IdxTrials6 = [];
    
    pos = -1;

    EV = 3;
    for i = 1:numel(Unit(NeuronIdx).Cnd)
        cnd = Unit(NeuronIdx).Cnd(i);
        if cnd.Value == [EV, pos]
            IdxTrials1 = [IdxTrials1; cnd.TrialIdx];
        end
    end

    EV = 6;
    for i = 1:numel(Unit(NeuronIdx).Cnd)
        cnd = Unit(NeuronIdx).Cnd(i);
        if cnd.Value == [EV, pos]
            IdxTrials2 = [IdxTrials2; cnd.TrialIdx];
        end
    end

    EV = 9;
    for i = 1:numel(Unit(NeuronIdx).Cnd)
        cnd = Unit(NeuronIdx).Cnd(i);
        if cnd.Value == [EV, pos]
            IdxTrials3 = [IdxTrials3; cnd.TrialIdx];
        end
    end
    
    pos = 1;
    
    EV = 3;
    for i = 1:numel(Unit(NeuronIdx).Cnd)
        cnd = Unit(NeuronIdx).Cnd(i);
        if cnd.Value == [EV, pos]
            IdxTrials4 = [IdxTrials4; cnd.TrialIdx];
        end
    end

    EV = 6;
    for i = 1:numel(Unit(NeuronIdx).Cnd)
        cnd = Unit(NeuronIdx).Cnd(i);
        if cnd.Value == [EV, pos]
            IdxTrials5 = [IdxTrials5; cnd.TrialIdx];
        end
    end

    EV = 9;
    for i = 1:numel(Unit(NeuronIdx).Cnd)
        cnd = Unit(NeuronIdx).Cnd(i);
        if cnd.Value == [EV, pos]
            IdxTrials6 = [IdxTrials6; cnd.TrialIdx];
        end
    end

    y = [];
    CntAll = [];
    for i = 1:numel(Unit(NeuronIdx).Trls)
        data = Unit(NeuronIdx).Trls(i);
        [cnts,~] = PSTH(data, WidthWindow, NumBins, cntrs);
        CntAll = [CntAll; cnts];
        if ~isempty(find(IdxTrials1==i))
            y = [y; 1];
        elseif ~isempty(find(IdxTrials2==i))
            y = [y; 2];
        elseif ~isempty(find(IdxTrials3==i))
            y = [y; 3];
        elseif ~isempty(find(IdxTrials4==i))
            y = [y; 4];
        elseif ~isempty(find(IdxTrials5==i))
            y = [y; 5];
        elseif ~isempty(find(IdxTrials6==i))
            y = [y; 6];
        end
    end
    
    if ShuffleStatus
        y = y(VecIdx);
    end
    
    mdl = fitglm(CntAll, y);
    pVal = coefTest(mdl);
    pValues = [pValues, pVal];
end

[~, col] = find(pValues<0.01);
BestUnits = NeuronIdxs(col);

figure('Units', 'normalized', 'OuterPosition', [0 0 1 1])
pltIdx = 1;
for NeuronIdx = BestUnits(randi(length(BestUnits), 12, 1))
    subplot(3, 4, pltIdx)
    plot_PSTH(Unit, NeuronIdx, WidthWindow, NumBins, cntrs)
    legend('off')
    title("PSTH of Unit " + num2str(NeuronIdx) + " (pValue = " + num2str(pValues(NeuronIdx) + ")", 4))
    pltIdx = pltIdx+1;
end

legend({'[3 -1]', '[3 1]', '[6 -1]', '[6 1]', '[9 -1]', '[9 1]', 'Average'})
REV_AllCndpValuesMean = mean(pValues)

set(gcf, 'PaperPositionMode', 'auto')
print("Step(2_3)_All_1_ShuffleStatus=" + num2str(ShuffleStatus), '-dpng', '-r0')

figure
histogram(pValues, 'Normalization', 'pdf')
title('Probability Density Function of p-Values for Cue Position')
xlabel('pValue')
ylabel('Probability')
xlim([0 1])

set(gcf, 'PaperPositionMode', 'auto')
print("Step(2_3)_All_2_ShuffleStatus=" + num2str(ShuffleStatus), '-dpng', '-r0')

figure
hold on
CntAllMean = [];
for CueValue = 3:3:9
    for pos = [-1, 1]
        CntAll = [];
        value = [CueValue, pos];
        
        for NeuronIdx = BestUnits
            data = GetCnd(Unit, NeuronIdx, value);
            [cnts,~] = PSTH(data, WidthWindow, NumBins, cntrs);
            CntAll = [CntAll; cnts];
        end
        
        CntAll = mean(CntAll, 1);
        CntAllMean = [CntAllMean; CntAll];
        plot(cntrs, CntAll, 'LineWidth', 1)
    end
end
CntAllMean = mean(CntAllMean, 1);

plot(cntrs, CntAllMean, 'k', 'LineWidth', 1.5)
title('Average PSTH of The Units (pValue < 0.01)');
xline(0, '--', 'Reward Cue', LabelHorizontalAlignment = 'center');
xline(0.3, '--', 'Delay Period', LabelHorizontalAlignment = 'center');
xline(0.9, '--', 'Reaction', LabelHorizontalAlignment = 'center');
xlim([-1.2, 2])
xlabel("Time (sec)")
ylabel('Firing Rate (Hz)')
legend('[3 -1]', '[3 1]', '[6 -1]', '[6 1]', '[9 -1]', '[9 1]', 'Avg')

set(gcf, 'PaperPositionMode', 'auto')
print("Step(2_3)_All_3_ShuffleStatus=" + num2str(ShuffleStatus), '-dpng', '-r0')
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