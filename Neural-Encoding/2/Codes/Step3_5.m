clc;clear;close all;

%% Loading data and initializing
load 'UnitsData.mat';
NumBins = 64; % number of bins
WidthWindow = 9; % width of the window
cntrs = linspace(-1.2, 2, NumBins);
CntAll = [];
AllDataTmp = zeros(numel(Unit), NumBins, 6);
AllData = zeros(2*numel(Unit), NumBins, 3);
AllDataReduced = zeros(2, size(AllData, 2), size(AllData, 3));
CueIdx = 1;

%% Step 3_5 - Just Expected Value - Reducing units dimensions to 2
for CueValue = 3:3:9
    for pos = [-1, 1]
        value = [CueValue, pos];
        for NeuronIdx = 1:numel(Unit)
            data = GetCnd(Unit, NeuronIdx, value);
            [cnts,~] = PSTH(data, WidthWindow, NumBins, cntrs);
            AllDataTmp(NeuronIdx, :, CueIdx) = cnts;
        end
        CueIdx = CueIdx+1;
    end
    AllData(:,:,CueValue/3) = [AllDataTmp(:,:,CueIdx-2); AllDataTmp(:,:,CueIdx-1)];
end

figure
for i = 1:3
    cov_mat = cov(AllData(:,:,i)');
    [V,D] = eig(cov_mat);
    D = diag(D);
    [D, I] = sort(D, 'descend');
    D = diag(D);
    V = V(:, I);
    V(:,3:end) = 0;

    B = sqrt(inv(D));
    A = V';
    Z = B*A*AllData(:,:,i);

    AllDataReduced(:,:,i) = Z(1:2, :);
    plot3(cntrs, AllDataReduced(1,:,i), AllDataReduced(2,:,i), 'LineWidth', 1.5)
    hold on
end

title("Dimension-Reduced Activity's PSTH")
xlabel('Time (sec)')
ylabel('Dim-1')
zlabel('Dim-2')
xlim([-1.2, 2])
[Y,Z] = meshgrid(min(AllDataReduced(1,:,:),[],'all'):0.1:max(AllDataReduced(1,:,:),[],'all'),min(AllDataReduced(2,:,:),[],'all'):0.1:max(AllDataReduced(2,:,:),[],'all'));
surf(Z*0,Y,Z, 'FaceAlpha', 0.1, 'FaceColor', 'b', 'EdgeColor', 'none');
surf(Z*0+0.3,Y,Z, 'FaceAlpha', 0.1, 'FaceColor', 'r', 'EdgeColor', 'none');
surf(Z*0+0.9,Y,Z, 'FaceAlpha', 0.1, 'FaceColor', 'k', 'EdgeColor', 'none');
legend('3', '6', '9', 'Reward Cue', 'Delay Period', 'Reaction')
hold off
set(gcf, 'PaperPositionMode', 'auto')
print("Step(3_5)", '-dpng', '-r0')

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