clc; clear; close all;

%% At first, run the 'script_mean_firing_rates_gratings_NEW.m' file
% In this section, I have removed Step 2 from the original file to include all neurons.

%% Load the 'S_monkey{1,2,3}.mat' files from the 'spikes_gratings' folder
for num_monkey = 1:3
    S_data{num_monkey} = load("./pvc-11/data_and_scripts/spikes_gratings/S_monkey" + num2str(num_monkey) + ".mat");
    monkey_data{num_monkey} = load("./pvc-11/data_and_scripts/spikes_gratings/data_monkey" + num2str(num_monkey) + "_gratings.mat");
end

%% Q1 - PSTHs
maxActivity = zeros(size(S_data, 2), size(S_data{1}.S, 2));
for num_monkey = 1:size(S_data, 2)
    for num_grating = 1:size(S_data{1}.S, 2)
        act = mean(S_data{num_monkey}.S(num_grating).mean_FRs, 2);
        maxActivity(num_monkey, num_grating) = find(act == max(act));
    end
end


for num_grating = 1:size(S_data{1}.S, 2)
    figure
    for num_monkey = 1:size(S_data, 2)
        c = [0, 0, 0];
        c(num_monkey) = 1;
        activity = S_data{num_monkey}.S(num_grating).mean_FRs/20e-3; % each bin is 20 ms

        hold on
        plt1 = plot(0:20:980, activity(maxActivity(num_monkey, num_grating), :), 'color', c,  'LineWidth', 1.5);
        plt2 = plot(0:20:980, mean(activity, 1), 'color', [c .5]);
    end
    title("PSTH , " + num2str((num_grating-1)*30) + "°")
    legend("Neuron " + num2str(maxActivity(1, num_grating)) + " , Monkey 1", "Average , Monkey 1", ...
           "Neuron " + num2str(maxActivity(2, num_grating)) + " , Monkey 2", "Average , Monkey 2", ...
           "Neuron " + num2str(maxActivity(3, num_grating)) + " , Monkey 3", "Average , Monkey 3")
    xlabel('Time (ms)')
    ylabel('Firing Rate (Hz)')
    xlim([0, 980])
    ylim([0, 140])
    set(gcf, 'PaperPositionMode', 'auto')
    print("Results/Q1_PSTH_" + num2str((num_grating-1)*30) + "°", '-dpng', '-r0')
end

%% Q1 - Tuning Curves
clc; close all;
for num_monkey = 1:size(S_data, 2)
    neurons = unique(maxActivity(num_monkey, :));
    activity = zeros(size(neurons, 2), size(S_data{1}.S, 2));
    leg = strings(size(neurons, 2));
    figure
    for num_neuron = 1:size(neurons, 2)
        neuron = neurons(num_neuron);
        leg(num_neuron) = leg(num_neuron) + "Neuron " + num2str(neuron);
        for num_grating = 1:size(S_data{1}.S, 2)
            activity(num_neuron, num_grating) = mean(S_data{num_monkey}.S(num_grating).mean_FRs(neuron,:))/20e-3; % each bin is 20 ms
        end
    end
    hold on
    plt = plot(0:30:330, activity, 'LineWidth', 1.5);
    title("Tuning Curves, Monkey " + num2str(num_monkey))
    xlabel('Orientation')
    ylabel('Firing Rate (Hz)')
    xlim([0, 330])
    ylim([0, 70])
    legend(leg)
    set(gcf, 'PaperPositionMode', 'auto')
    print("Results/Q1_TC_Monkey" + num2str(num_monkey), '-dpng', '-r0')
end

%% Q2 - Preferred Orientations, (0° - 330°)
clc; close all;
for num_monkey = 1:size(S_data, 2)
    channels = monkey_data{num_monkey}.data.CHANNELS;
    map = monkey_data{num_monkey}.data.MAP;
    neurons = 1:size(S_data{num_monkey}.S(1).mean_FRs, 1);
    PrfGrating = zeros(1, length(neurons));
    PrfGratingTmp = zeros(12, length(neurons));
    
    for num_grating = 1:size(S_data{1}.S, 2)
        act = mean(S_data{num_monkey}.S(num_grating).mean_FRs, 2)';
        PrfGratingTmp(num_grating, :) = act;
    end
    
    for neuron = neurons
        temp = find(PrfGratingTmp(:, neuron) == max(PrfGratingTmp(:, neuron)));
        PrfGrating(neuron) = temp(1);
    end
    
    GratingMap = nan(size(map));
    for i = 1:size(map, 1)
        for j = 1:size(map, 2)
            if ~isnan(map(i, j))
                temp = find(channels(:, 1) == map(i, j));
                if ~isempty(temp)
                    GratingMap(i, j) = (PrfGrating(temp(1))-1)*30;
                end
            end
        end
    end
    FinalGrating{num_monkey} = GratingMap;
    
    % plot
    figure
    % FinalGrating{num_monkey}(isnan(FinalGrating{num_monkey})) = 0;
    plt = imagesc(FinalGrating{num_monkey});
    title("Preferred Orientations, (0° - 330°), " + "Monkey " + num2str(num_monkey))
    caxis([0, 330])
    colorbar
    set(plt, 'AlphaData', ~isnan(FinalGrating{num_monkey}))
    colormap jet
    
    set(gcf, 'PaperPositionMode', 'auto')
    print("Results/Q2_PO1(0°-330°)_Monkey" + num2str(num_monkey), '-dpng','-r0')
    
    figure
    FinalGrating{num_monkey}(isnan(FinalGrating{num_monkey})) = 0;
    pcolor(1:10, 1:10, FinalGrating{num_monkey});
    title("Preferred Orientations, (0° - 330°), " + "Monkey " + num2str(num_monkey))
    caxis([0, 330])
    colorbar
    shading interp
    axis ij
    colormap jet

    set(gcf, 'PaperPositionMode', 'auto')
    print("Results/Q2_PO2(0°-330°)_Monkey" + num2str(num_monkey), '-dpng', '-r0')
end

%% Q2 - Preferred Orientations, (0° - 150°)
clc; close all;
for num_monkey = 1:size(S_data, 2)
    channels = monkey_data{num_monkey}.data.CHANNELS;
    map = monkey_data{num_monkey}.data.MAP;
    neurons = 1:size(S_data{num_monkey}.S(1).mean_FRs, 1);
    PrfGrating = zeros(1, length(neurons));
    PrfGratingTmp = zeros(12, length(neurons));
    
    for num_grating = 1:size(S_data{1}.S, 2)
        act = mean(S_data{num_monkey}.S(num_grating).mean_FRs, 2)';
        PrfGratingTmp(num_grating, :) = act;
    end
    
    for neuron = neurons
        temp = find(PrfGratingTmp(:, neuron) == max(PrfGratingTmp(:, neuron)));
        PrfGrating(neuron) = temp(1);
    end
    
    GratingMap = nan(size(map));
    for i = 1:size(map, 1)
        for j = 1:size(map, 2)
            if ~isnan(map(i, j))
                temp = find(channels(:, 1) == map(i, j));
                if ~isempty(temp)
                    GratingMap(i, j) = (PrfGrating(temp(1))-1)*30;
                end
            end
            % 0°-150°
            if GratingMap(i, j) > 150
                GratingMap(i, j) = GratingMap(i, j) - 180;
            end
        end
    end
    FinalGrating{num_monkey} = GratingMap;

    % plot
    figure
    % FinalGrating{num_monkey}(isnan(FinalGrating{num_monkey})) = 0;
    plt = imagesc(FinalGrating{num_monkey});
    title("Preferred Orientations, (0° - 150°), " + "Monkey " + num2str(num_monkey))
    caxis([0, 150])
    colorbar
    set(plt, 'AlphaData', ~isnan(FinalGrating{num_monkey}))
    colormap jet
    
    set(gcf, 'PaperPositionMode', 'auto')
    print("Results/Q2_PO3(0°-150°)_Monkey" + num2str(num_monkey), '-dpng','-r0')
    
    figure
    FinalGrating{num_monkey}(isnan(FinalGrating{num_monkey})) = 0;
    pcolor(1:10, 1:10, FinalGrating{num_monkey});
    title("Preferred Orientations, (0° - 150°), " + "Monkey " + num2str(num_monkey))
    caxis([0, 150])
    colorbar
    shading interp
    axis ij
    colormap jet

    set(gcf, 'PaperPositionMode', 'auto')
    print("Results/Q2_PO4(0°-150°)_Monkey" + num2str(num_monkey), '-dpng', '-r0')
end

%% Q3 - find the dependence of r_sc (noise correlation)
clc; close all;
ShuffleStatus = input('Set status of the shuffling (Q3): (1=on, 0=off) '); % shuffle the final labels
for num_monkey = 1:size(S_data, 2)
    neurons = 1:size(S_data{num_monkey}.S(1).mean_FRs, 1);
    activity = zeros(size(neurons, 2), size(S_data{1}.S, 2));
    
    % calculate tuning curves
    for num_neuron = 1:size(neurons, 2)
        neuron = neurons(num_neuron);
        for num_grating = 1:size(S_data{1}.S, 2)
            activity(num_neuron, num_grating) = mean(S_data{num_monkey}.S(num_grating).mean_FRs(neuron,:))/20e-3; % each bin is 20 ms
        end
    end
    TuningCurves{num_monkey} = activity;
    
    % calculate r_s
    r_s_tmp = zeros(size(activity, 1));
    for n1 = 1:size(activity, 1)-1
        for n2 = n1+1:size(activity, 1)
            corr_tmp = corrcoef(activity(n1, :), activity(n2, :));
            r_s_tmp(n1, n2) = corr_tmp(2,1);
        end
    end
    r_s{num_monkey} = r_s_tmp;
    
    % calculate r_sc
    S_data_copy = S_data;
    firing_rates = zeros(size(S_data_copy{1}.S, 2) * size(S_data_copy{num_monkey}.S(1).trial, 2), size(neurons, 2));
    for num_neuron = 1:size(neurons, 2)
        fr_tmp = zeros(size(S_data_copy{1}.S, 2) * size(S_data_copy{num_monkey}.S(1).trial, 2), 1);
        for num_grating = 1:size(S_data_copy{1}.S, 2)
            if ShuffleStatus
                S_data_copy{num_monkey}.S(1).trial = S_data_copy{num_monkey}.S(1).trial(randperm(size(S_data_copy{num_monkey}.S(1).trial, 2)));
            end

            for num_trial = 1:size(S_data_copy{num_monkey}.S(1).trial, 2)
                fr_tmp((num_grating-1) * size(S_data_copy{num_monkey}.S(1).trial, 2) + num_trial, 1) = sum(S_data_copy{num_monkey}.S(1).trial(num_trial).counts(num_neuron, :));
            end
            fr_tmp((num_grating-1) * size(S_data_copy{num_monkey}.S(1).trial, 2) + 1 :  num_grating * size(S_data_copy{num_monkey}.S(1).trial, 2)) = zscore(fr_tmp((num_grating-1) * size(S_data{num_monkey}.S(1).trial, 2) + 1 :  num_grating * size(S_data{num_monkey}.S(1).trial, 2)));
        end
        firing_rates(:, num_neuron) = fr_tmp;
    end
    r_sc_tmp = corrcoef(firing_rates);
    r_sc{num_monkey} = r_sc_tmp;
    
    % calculate distance between neurons
    distance_channels = 400e-6; % The distance on the array between adjacent channels is 400um.
    distance_tmp = zeros(size(S_data{1}.S(1).trial(1).spikes, 1));
    channels = monkey_data{num_monkey}.data.CHANNELS;
    map = monkey_data{num_monkey}.data.MAP;  
    for n1 = 1:size(neurons, 2)-1
        for n2 = n1+1:size(neurons, 2)
            inf1 = channels(n1, 1);
            inf2 = channels(n2, 1);
            [r1, c1] = find(inf1 == map);
            [r2, c2] = find(inf2 == map);
            dx = abs(r1 - r2);
            dy = abs(c1 - c2);
            distance_tmp(n1, n2) = distance_channels * sqrt(dx^2 + dy^2);
        end
    end
    distance{num_monkey} = distance_tmp;

    % plot: distance = f(r_sc) 
    r_s_vec = [];
    r_sc_vec = [];
    distance_vec = [];
    for n1 = 1:size(neurons, 2)-1
        for n2 = n1+1:size(neurons, 2)
            r_s_vec = [r_s_vec, r_s_tmp(n1, n2)];
            r_sc_vec = [r_sc_vec, r_sc_tmp(n1, n2)];
            distance_vec = [distance_vec, distance_tmp(n1, n2)];
        end
    end
    [distance_vec, ind] = sort(distance_vec);
    r_s_vec = r_s_vec(ind);
    r_sc_vec = r_sc_vec(ind);
    
    Matrix = zeros(4, 17, 3);
    cnt = 1;
    for dis = 0.25e-3:0.25e-3:4.25e-3
        [~, col1] = find(distance_vec >= dis-0.26e-3 & distance_vec < dis+0.26e-3);
        r_s_value = r_s_vec(col1);
        r_sc_value = r_sc_vec(col1);
        distance_value = distance_vec(col1);
       
        [~, col2] = find(r_s_value >= 0.5);
        r_sc_value2 = r_sc_value(col2);
        r_sc_mean = mean(r_sc_value2);
        r_sc_var = var(r_sc_value2);
        distance_value2 = distance_value(col2);
        distance_mean = mean(distance_value2);
        distance_var = var(distance_value2);
        Matrix(1, cnt, :) = [distance_mean, r_sc_mean, r_sc_var];

        [~, col2] = find(r_s_value >= 0 & r_s_value < 0.5);
        r_sc_value2 = r_sc_value(col2);
        r_sc_mean = mean(r_sc_value2);
        r_sc_var = var(r_sc_value2);
        distance_value2 = distance_value(col2);
        distance_mean = mean(distance_value2);
        distance_var = var(distance_value2);
        Matrix(2, cnt, :) = [distance_mean, r_sc_mean, r_sc_var];

        [~, col2] = find(r_s_value >= -0.5 & r_s_value < 0);
        r_sc_value2 = r_sc_value(col2);
        r_sc_mean = mean(r_sc_value2);
        r_sc_var = var(r_sc_value2);
        distance_value2 = distance_value(col2);
        distance_mean = mean(distance_value2);
        distance_var = var(distance_value2);
        Matrix(3, cnt, :) = [distance_mean, r_sc_mean, r_sc_var];

        [~, col2] = find(r_s_value < -0.5);
        r_sc_value2 = r_sc_value(col2);
        r_sc_mean = mean(r_sc_value2);
        r_sc_var = var(r_sc_value2);
        distance_value2 = distance_value(col2);
        distance_mean = mean(distance_value2);
        distance_var = var(distance_value2);
        Matrix(4, cnt, :) = [distance_mean, r_sc_mean, r_sc_var];

        cnt = cnt + 1;
    end

    figure
    hold on
    for i = 1:4
        plot(1000 * Matrix(i, :, 1), Matrix(i, :, 2), 'LineWidth', 5 - i, 'color', [0, 0, 0] + i/6);
    end
    for i = 1:4
        errorbar(1000 * Matrix(i, :, 1), Matrix(i, :, 2), Matrix(i, :, 3), 'color', [0, 0, 0] + i/6);
    end
    title("Monkey " + num2str(num_monkey))
    xlabel("Distance (mm)")
    ylabel("r_{sc}")
    xlim([0, 4.5])
    ylim([-0.01, 0.4])
    legend('0.5 <= r_s', '0 <= r_s < 0.5', '-0.5 <= r_s < 0', 'r_s < -0.5', 'NumColumns', 1)
    set(gcf, 'PaperPositionMode', 'auto')
    print("Results/Q3_1_(rsc_distance)_Monkey" + num2str(num_monkey) + "_Shuffle=" + num2str(ShuffleStatus), '-dpng', '-r0')

% plot: r_sc = f(r_s) 
    r_s_vec = [];
    r_sc_vec = [];
    distance_vec = [];
    for n1 = 1:size(neurons, 2)-1
        for n2 = n1+1:size(neurons, 2)
            r_s_vec = [r_s_vec, r_s_tmp(n1, n2)];
            r_sc_vec = [r_sc_vec, r_sc_tmp(n1, n2)];
            distance_vec = [distance_vec, distance_tmp(n1, n2)];
        end
    end
    [r_s_vec, ind] = sort(r_s_vec);
    r_sc_vec = r_sc_vec(ind);
    distance_vec = distance_vec(ind);
    
    Matrix = zeros(4, 21, 3);
    cnt = 1;
    for rs = -1:0.1:1
        [~, col1] = find(r_s_vec >= rs-0.126 & r_s_vec < rs+0.126);
        r_s_value = r_s_vec(col1);
        r_sc_value = r_sc_vec(col1);
        distance_value = distance_vec(col1);
       
        [~, col2] = find(distance_value >= 3e-3);
        r_sc_value2 = r_sc_value(col2);
        r_sc_mean = mean(r_sc_value2);
        r_sc_var = var(r_sc_value2);
        r_s_value2 = r_s_value(col2);
        r_s_mean = mean(r_s_value2);
        r_s_var = var(r_s_value2);
        Matrix(1, cnt, :) = [r_s_mean, r_sc_mean, r_sc_var];

        [~, col2] = find(distance_value >= 2e-3 & distance_value < 3e-3);
        r_sc_value2 = r_sc_value(col2);
        r_sc_mean = mean(r_sc_value2);
        r_sc_var = var(r_sc_value2);
        r_s_value2 = r_s_value(col2);
        r_s_mean = mean(r_s_value2);
        r_s_var = var(r_s_value2);
        Matrix(2, cnt, :) = [r_s_mean, r_sc_mean, r_sc_var];

        [~, col2] = find(distance_value >= 1e-3 & distance_value < 2e-3);
        r_sc_value2 = r_sc_value(col2);
        r_sc_mean = mean(r_sc_value2);
        r_sc_var = var(r_sc_value2);
        r_s_value2 = r_s_value(col2);
        r_s_mean = mean(r_s_value2);
        r_s_var = var(r_s_value2);
        Matrix(3, cnt, :) = [r_s_mean, r_sc_mean, r_sc_var];

        [~, col2] = find(distance_value < 1e-3 & distance_value >= 0);
        r_sc_value2 = r_sc_value(col2);
        r_sc_mean = mean(r_sc_value2);
        r_sc_var = var(r_sc_value2);
        r_s_value2 = r_s_value(col2);
        r_s_mean = mean(r_s_value2);
        r_s_var = var(r_s_value2);
        Matrix(4, cnt, :) = [r_s_mean, r_sc_mean, r_sc_var];

        cnt = cnt + 1;
    end

    figure
    hold on
    for i = 1:4
        plot(Matrix(i, :, 1), Matrix(i, :, 2), 'LineWidth', 5 - i, 'color', [0, 0, 0] + i/6);
    end
    for i = 1:4
        errorbar(Matrix(i, :, 1), Matrix(i, :, 2), Matrix(i, :, 3), 'color', [0, 0, 0] + i/6);
    end
    title("Monkey " + num2str(num_monkey))
    xlabel("r_s")
    ylabel("r_{sc}")
    xlim([-1, 1])
    ylim([-0.05, 0.45])
    legend('3mm <= Distance', '2mm <= Distance < 3mm', '1mm <= Distance < 2mm', '0 <= Distance < 1mm', 'NumColumns', 1)
    set(gcf, 'PaperPositionMode', 'auto')
    print("Results/Q3_2_(rsc_rs)_Monkey" + num2str(num_monkey) + "_Shuffle=" + num2str(ShuffleStatus), '-dpng', '-r0')
    
    % plot: r_s = f(distance)
    rs = -1:0.1:1;
    dis = 0.25e-3:0.25e-3:4.25e-3;
    r_s_vec = [];
    r_sc_vec = [];
    distance_vec = [];
    for n1 = 1:size(neurons, 2)-1
        for n2 = n1+1:size(neurons, 2)
            r_s_vec = [r_s_vec, r_s_tmp(n1, n2)];
            r_sc_vec = [r_sc_vec, r_sc_tmp(n1, n2)];
            distance_vec = [distance_vec, distance_tmp(n1, n2)];
        end
    end

    Matrix = zeros(length(rs), length(dis));
    Values = [r_s_vec; distance_vec; r_sc_vec];
    cnt_dis = 1;
    for d = dis
        cnt_rs = 1;
        for r = rs
            [~, col] = find(Values(1, :) > r-0.5 & Values(1, :) < r+0.5 & Values(2, :) > d-0.125e-3 & Values(2, :) < d+0.125e-3);
            if ~isempty(col)
                Matrix(cnt_rs, cnt_dis) = Values(3, col(1));
            end
            cnt_rs = cnt_rs + 1;
        end
        cnt_dis = cnt_dis + 1;
    end
    figure
    pcolor(1000 * dis, rs, Matrix);
    title("Monkey " + num2str(num_monkey));
    xlabel("Distance (mm)");
    ylabel("r_s");
    a = colorbar;
    ylabel(a, 'r_{sc}', 'Rotation', 0, 'Position', [0.5 -0.05]);
    caxis([-0.05, 0.45])
    colormap jet
    shading interp
    set(gcf, 'PaperPositionMode', 'auto')
    print("Results/Q3_3_(rs_rsc_distance)_Monkey" + num2str(num_monkey) + "_Shuffle=" + num2str(ShuffleStatus), '-dpng', '-r0')
end