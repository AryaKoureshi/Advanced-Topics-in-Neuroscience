clc; clear; close all;

%% =========== part 1 ===========
%% Q2
numExperiment = 20;
sigma = 1;
dt = 0.1;
timeInterval = 0:dt:10;
biases = [0 0.1 1 10];

for bb = 1:size(biases, 2)
    bias = biases(bb);
    X = zeros(length(timeInterval), numExperiment);
    choices = zeros(1, numExperiment);
    
    for iter = 1:numExperiment
        [X(:, iter), choices(1, iter)] = simple_model(bias, sigma, dt, timeInterval);
    end
    
    randomData = normrnd(bias*timeInterval(end), sigma*timeInterval(end), size(X,1), size(X,2));
    figure
    histogram(X, 'Normalization', 'pdf', 'FaceColor','red')
    hold on 
    histogram(randomData, 'Normalization', 'pdf', 'FaceColor', 'blue')
    xlabel('X')
    ylabel('Probability Density')
    title("Distribution, Bias = " + num2str(bias))
    legend('X', 'N(0, \sigma)')
    set(gcf, 'PaperPositionMode', 'auto')
    print("Results/pdfBias" + num2str(bias),'-dpng','-r0')

    figure
    plot(timeInterval, mean(X, 2), 'k', 'LineWidth', 2)
    hold on
    plot(timeInterval, X', '--')
    legend('Average of all Experiments', 'Different Experiments')
    title("Bias = " + num2str(bias))
    xlabel('Time (sec)')
    ylabel('X')
    set(gcf, 'PaperPositionMode', 'auto')
    print("Results/different_Xs_bias_" + num2str(bias),'-dpng','-r0')
end

X = zeros(length(biases), length(timeInterval));
for bb = 1:size(biases, 2)
    bias = biases(bb);
    [X(bb, :), ~] = simple_model(bias, sigma, dt, timeInterval);
end

figure
hold on
for i = 1:size(X, 1)
    plot(timeInterval, X(i,:))
end
xlabel('Time (sec)')
ylabel('X')
legend('B = 0', 'B = 0.1', 'B = 1', 'B = 10')
set(gcf, 'PaperPositionMode', 'auto')
print("Results/10sec_different_biases",'-dpng','-r0')

%% Q3
numExperiment = 1000;
sigma = 1;
dt = 0.1;
bias = 0.1;
times = 0.5:1:50;
error = zeros(1, length(times));
errorTheory = zeros(1, length(times));

cnt = 1;
for time = times
    timeInterval = 0:dt:time;
    choices = zeros(1, numExperiment);
    X = zeros(length(timeInterval), numExperiment);
    for iter = 1:numExperiment
        [X(:, iter), choices(1, iter)] = simple_model(bias, sigma, dt, timeInterval);
    end

    error(1, cnt) = 1-sum(choices == sign(bias))/length(choices);
    meanTheory = bias*time;
    sigmaTheory = sqrt(dt*time);
    errorProb = 1 - min(8*sigmaTheory, meanTheory + 4*sigmaTheory)/(8*sigmaTheory);
    errorTheory(1, cnt) = errorProb;
    cnt = cnt + 1;
end

plot(times, error)
hold on
plot(times, errorTheory)
xlabel('Time Intervals')
ylabel('Error')
legend('Simulation', 'Theory')
ylim([0, 0.6])

set(gcf, 'PaperPositionMode', 'auto')
print("Results/Error",'-dpng','-r0')

close all;
x = -5:.1:5;
y = normpdf(x,0,1);
plot(x, y)
xl = xline(0, '--', '\mu');
xl = xline(-4, '--', '\mu + -4\sigma');
xl = xline(4, '--', '\mu + 4\sigma');
xlim([-5 5])
ylim([0, 0.5])
set(gca,'xtick',[])
set(gca,'ytick',[])

set(gcf, 'PaperPositionMode', 'auto')
print('Results/normalProof','-dpng','-r0')

%% Q4
numExperiment = 10;
sigma = 1;
dt = 0.1;
timeInterval = 0:dt:10;
bias = 0.1;

X = zeros(length(timeInterval), numExperiment);
for iter = 1:numExperiment
    [X(:, iter), ~] = simple_model(bias, sigma, dt, timeInterval);
end

mean_X = mean(X, 2);
var_X = var(X, 1, 2);

subplot(1,2,1)
plot(timeInterval, X)
xlabel('Time (sec)')
title('X(t)')

subplot(1,2,2)
plot(timeInterval, mean_X, 'b')
hold on
plot(timeInterval, bias*timeInterval, '--b')
plot(timeInterval, var_X, 'r')
plot(timeInterval, sigma*timeInterval, '--r')
legend('Simulation (Mean)', 'Theory (Mean)', 'Simulation (Variance)', 'Theory (Variance)')
xlabel('Time (sec)')
title('Mean and Variance')

set(gcf, 'PaperPositionMode', 'auto')
print("Results/X_MeanVar",'-dpng','-r0')

%% Q5
bias = 0;
sigma = 1;
Xs0 = -10:1:10;
timeLimits = 1:1:100;

probs = zeros(length(timeLimits), length(Xs0));

cnt = 1;
for X0 = Xs0
    cnt2 = 1;
    for timeLimit = timeLimits
        p = simple_model2(bias, sigma, X0, timeLimit);
        probs(cnt2, cnt) = p;
        cnt2 = cnt2+1;
    end
    cnt = cnt+1;
end

figure
hold on
for cnt2 = 1:length(timeLimits)
    plot(Xs0, probs(cnt2, :))
end
xlabel('Start Point')
ylabel('Probability')
ylim([0, 1])
colormap default
c = colorbar('Ticks', [timeLimits(1), timeLimits(end/2), timeLimits(end)], 'TickLabels', ...
        {num2str(timeLimits(1)), num2str(timeLimits(end/2)), num2str(timeLimits(end))});
c.Label.String = 'Time Limit (sec)';
caxis([timeLimits(1), timeLimits(end)])

set(gcf, 'PaperPositionMode', 'auto')
print("Results/Q05_plot",'-dpng','-r0')

figure
imagesc(Xs0, timeLimits, probs)
xlabel('Start Point')
ylabel('Time Limit (sec)')
colormap default
c = colorbar;
c.Label.String = 'Probability';

set(gcf, 'PaperPositionMode', 'auto')
print("Results/Q05_countorPlot",'-dpng','-r0')

%% Q7
numExperiment = 1000;
thresholds = [-10, 10];
bias = 0.1;
sigma = 1;
X0 = 0;
dt = 0.01;

ts = zeros(1, numExperiment);
choices = zeros(1, numExperiment);
for iter_no = 1:numExperiment
    [t, choice] = two_choice_trial(thresholds, bias, sigma, X0, dt);
    ts(1, iter_no) = t; 
    choices(1, iter_no) = choice; 
end

figure
histAll = histogram(ts, 'Normalization', 'pdf');
histAllBins = (histAll.BinEdges(2:end)+histAll.BinEdges(1:end-1))/2;
histAllValues = histAll.Values;

histCorrect = histogram(ts(choices==1), 'Normalization', 'pdf');
histCorrectBins = (histCorrect.BinEdges(2:end)+histCorrect.BinEdges(1:end-1))/2;
histCorrectValues = histCorrect.Values;

histWrong = histogram(ts(choices==-1), 'Normalization', 'pdf');
histWrongBins = (histWrong.BinEdges(2:end)+histWrong.BinEdges(1:end-1))/2;
histWrongValues = histWrong.Values;

plot(histAllBins/dt, histAllValues*dt)
hold on
pd = makedist('InverseGaussian', 'mu', thresholds(2)/bias, 'lambda', (thresholds(2)/sigma)^2);
pdf_values = pdf(pd, histAllBins/dt);
plot(histAllBins/dt, pdf_values)
xlabel('Reaction Time')
ylabel('Probability Density')
legend('Simulation', 'Theory')

set(gcf, 'PaperPositionMode', 'auto')
print("Results/Q07_All",'-dpng','-r0')

figure
hold on
plot(histCorrectBins/dt, histCorrectValues*dt)
plot(histWrongBins/dt, histWrongValues*dt)
xlabel('Reaction Time')
ylabel('Probability Density')
legend('Correct choices', 'Wrong choices')

set(gcf, 'PaperPositionMode', 'auto')
print("Results/Q07_CorrectAndWrong",'-dpng','-r0')

%% Q8
threshold = 10;
biases = 1:1:10;
sigma = 1;
X_0 = [0;0];
dt = 0.01;
num_iters = 1000;
stats = zeros(length(biases), length(biases), num_iters);

bias_1_index = 1;
for bias_1 = biases
    bias_2_index = 1;
    for bias_2 = biases
        bias = [bias_1; bias_2];
        for iter_no = 1:num_iters
            stats(bias_1_index, bias_2_index, iter_no)  = race_trial(threshold, bias, sigma, X_0, dt);
        end
        bias_2_index = bias_2_index+1;
    end
    bias_1_index = bias_1_index+1;
end

imagesc(biases, biases, 100*sum(stats==1, 3)/size(stats, 3))
set(gca,'YDir','normal')
xlabel('Bias (Choice 2)')
ylabel('Bias (Choice 1)')
colormap default
c = colorbar;
c.Label.String = 'Win Rate (Choice 1)';

set(gcf, 'PaperPositionMode', 'auto')
print("Results/Q08",'-dpng','-r0')

%% Q9
threshold = 10;
biases = 1:1:10;
sigma = 1;
X_0 = [0;0];
dt = 0.01;
num_iters = 1000;
time_limit = 2;
stats = zeros(length(biases), length(biases), num_iters);

bias_1_index = 1;
for bias_1 = biases
    bias_2_index = 1;
    for bias_2 = biases
        bias = [bias_1; bias_2];
        for iter_no = 1:num_iters
            stats(bias_1_index, bias_2_index, iter_no)  = race_trial_timeLimit(threshold, bias, sigma, X_0, dt, time_limit);
        end
        bias_2_index = bias_2_index+1;
    end
    bias_1_index = bias_1_index+1;
end

imagesc(biases, biases, 100*sum(stats==3, 3)/size(stats, 3))
set(gca,'YDir','normal')
xlabel('Bias (Choice 2)')
ylabel('Bias (Choice 1)')
colormap default
c = colorbar;
c.Label.String = 'Percentage of trials when time exceeds the threshold';

set(gcf, 'PaperPositionMode', 'auto')
print("Results/Q09",'-dpng','-r0')

%% ========== part 2 =========
%% Q1
LIPthreshold = 50;
MT_pValues = [0.1; 0.05];
LIP_weights = [0.1; -0.2];
[rt, lip_eventTimes, MT, times] = lip_activity(MT_pValues, LIP_weights, LIPthreshold);
figure
hold on
scatter(lip_eventTimes, lip_eventTimes*0+3, 'k', 'filled')
scatter(times(MT(1,:)==1), MT(1, MT(1,:)==1)*0+2, 'b', 'filled')
scatter(times(MT(2,:)==1), MT(2, MT(2,:)==1)*0+1, 'r', 'filled')
legend('LIP Neuron', 'MT Neuron excitatory', 'MT Neuron inhibitory')
ylim([0, 4])
xlim([lip_eventTimes(1)-0.01, lip_eventTimes(end)+0.01])
xlabel('Time (sec)')

set(gcf, 'PaperPositionMode', 'auto')
print("Results/part2_Q1",'-dpng','-r0')

%% Q2

sigma = 0.4;
stimuli_values = 0:0.01:6;
Mt_neuron_1 = makedist('Normal','mu',2,'sigma',sigma);
Mt_neuron_2 = makedist('Normal','mu',4,'sigma',sigma);
figure
hold on
plot(stimuli_values, pdf(Mt_neuron_1, stimuli_values))
plot(stimuli_values, pdf(Mt_neuron_2, stimuli_values))
ylim([0, 1])
xlabel('Stimuli Value')
ylabel('Probability of Firing')
legend('MT Neuron 1', 'MT Neuron 2')

set(gcf, 'PaperPositionMode', 'auto')
print("Results/part2_Q2_1",'-dpng','-r0')

LIP_weights = [0.1 -0.1; -0.1 0.1];
numExperiment = 100;
[LIP_spikes, MT, times] = lip_activity_enhanced(LIP_weights, stimuli_values, Mt_neuron_1, Mt_neuron_2, numExperiment);

figure
hold on
scatter(times(LIP_spikes(1,:)==1), LIP_spikes(LIP_spikes(1,:)==1)*0+3, 'k', 'filled')
scatter(times(LIP_spikes(2,:)==1), LIP_spikes(LIP_spikes(2,:)==1)*0+2, 'b', 'filled')
scatter(times(MT(1,:)==1), MT(1, MT(1,:)==1)*0+1, 'g', 'filled')
scatter(times(MT(2,:)==1), MT(2, MT(2,:)==1)*0, 'r', 'filled')
legend('LIP1 = MT1-MT2', 'LIP2 = -MT1+MT2', 'MT Neuron 1', 'MT Neuron 2')
ylim([0, 4])
xlim([times(1)-0.01, times(end)+0.01])
xlabel('Time (sec)')

set(gcf, 'PaperPositionMode', 'auto')
print("Results/part2_Q2_2",'-dpng','-r0')

LIP_weights = [0.1 0.1; 0.1 -0.1];
numExperiment = 100;
[LIP_spikes, MT, times] = lip_activity_enhanced(LIP_weights, stimuli_values, Mt_neuron_1, Mt_neuron_2, numExperiment);

figure
hold on
scatter(times(LIP_spikes(1,:)==1), LIP_spikes(LIP_spikes(1,:)==1)*0+3, 'k', 'filled')
scatter(times(LIP_spikes(2,:)==1), LIP_spikes(LIP_spikes(2,:)==1)*0+2, 'b', 'filled')
scatter(times(MT(1,:)==1), MT(1, MT(1,:)==1)*0+1, 'g', 'filled')
scatter(times(MT(2,:)==1), MT(2, MT(2,:)==1)*0, 'r', 'filled')
legend('LIP1 = MT1+MT2', 'LIP2 = MT1-MT2', 'MT Neuron 1', 'MT Neuron 2')
ylim([0, 4])
xlim([times(1)-0.01, times(end)+0.01])
xlabel('Time (sec)')

set(gcf, 'PaperPositionMode', 'auto')
print("Results/part2_Q2_3",'-dpng','-r0')

%% Functions
function [X, choice] = simple_model(bias, sigma, dt, timeInterval)
    dW = normrnd(0, sqrt(dt), 1, length(timeInterval));
    X = zeros(1, length(timeInterval));
    for i = 1:1:length(timeInterval)-1
        dX = bias*dt + sigma*dW(i);
        X(i+1) = X(i) + dX;
    end
    choice = sign(X(end));
end

function p = simple_model2(bias, sigma, X0, timeLimit)
        mu = bias * timeLimit;
        sigma = sqrt(sigma * timeLimit);
        p = normcdf(X0, mu, sigma);
        p = 1 - p;
end

% Q6
function [t, choice] = two_choice_trial(thresholds, bias, sigma, X0, dt)
    X = X0;
    t = 0;
    while X > thresholds(1) && X < thresholds(2)
        dW = normrnd(0, sigma, 1, 1);
        dX = bias*dt + sigma*dW;
        X = X + dX;
        t = t + dt;
    end
    choice = sign(X);
end

function winner_no = race_trial(threshold, bias, sigma, X0, dt)
    X = X0;
    t = 0;
    while X < threshold 
        dW = normrnd(0,sigma,2,1);
        dX = bias*dt + sigma*dW;
        X = X + dX;
        t = t + dt;
    end
    [winner_no, ~] = find(X > threshold);
    if length(winner_no) > 1
        winner_no = 0;
    end
end

% Q9
function winner_no = race_trial_timeLimit(threshold, bias, sigma, X_0, dt, time_limit)
    X = X_0;
    t = 0;
    while X < threshold & t < time_limit
        dW = normrnd(0,sigma,2,1);
        dX = bias*dt+sigma*dW;
        X = X + dX;
        t = t+dt;
    end
    [winner_no, ~] = find(X > threshold);
    if length(winner_no) > 1
        winner_no = 0;
    elseif t > time_limit
        winner_no = 3;
    end
end

function [rt, lip_eventTimes, MT, times] = lip_activity(MT_pValues, LIP_weights, LIPthreshold)
    dt = 0.001;
    N = [0; 0];
    rate = 0.0;
    lip_eventTimes = [];
    t = 0;
    M = 100;
    MT = [];
    times = [];
    while rate < LIPthreshold
        times = [times, t];
        t = t + dt;
        dN = rand(2, 1) < MT_pValues;
        MT = [MT, dN];
        N = N + dN;
        p_LIP = LIP_weights'*N;
        LIP_event = rand(1,1) < p_LIP;
        if LIP_event
            lip_eventTimes = [lip_eventTimes, t];
        end
        if  length(lip_eventTimes)>M
            rate = M/(t-lip_eventTimes(end-M+1));
        end
    end
    rt = t;
end

function [LIP_spikes, MT, times] = lip_activity_enhanced(LIP_weights, stimuli_values, MT1_dist, MT2_dist, num_iters)
    dt = 0.001;
    N = [0; 0];
    LIP_spikes = [];
    t = 0;
    MT = [];
    times = [];
    for iter_no = 1:num_iters
        stimuli = stimuli_values(randi(length(stimuli_values)));
        times = [times, t];
        t = t + dt;
        firing_probs = [pdf(MT1_dist, stimuli); pdf(MT2_dist, stimuli)];
        dN = rand(2, 1) < firing_probs;
        N = N + dN;
        MT = [MT, dN];
        p_LIP = N'*LIP_weights;
        LIP_event = rand(2,1) < p_LIP';
        LIP_spikes = [LIP_spikes, LIP_event];
    end
end