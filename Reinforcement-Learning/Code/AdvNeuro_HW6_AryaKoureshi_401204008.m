clc; close all; clear;
%% Question 1
map = zeros(15, 15);
reward = [2, 12];
cat = [8, 4];
values = -1 * ones(15, 15);
values(sub2ind(size(values), [2, 8], [12, 4])) = [10, -10];
qmatrix = zeros(15, 15, 4);
runNum = 10000;
eps = 0.1;
lr = 0.5;
gamma = 0.8;
runSaver = zeros(2, 11, 1000);
selectedRuns = [5, 212, 666];
i = 0;
while i <= runNum
    i = i + 1;
    state = randi([1, 15], 1, 2);
    cnt = 0;
    while ~(isequal(state, reward) || isequal(state, cat))
        cnt = cnt + 1;
        if rand() < eps
            direction = randi(4);
        else
            [~, direction] = max(qmatrix(state(1), state(2), :));
        end

        if direction == 1
            state1 = [state(1), min(state(2) + 1, 15)];
        elseif direction == 2
            state1 = [state(1), max(state(2) - 1, 1)];
        elseif direction == 3
            state1 = [max(state(1) - 1, 1), state(2)];
        elseif direction == 4
            state1 = [min(state(1) + 1, 15), state(2)];
        end

        rewardM = values(state1(1), state1(2));
        qmatrix(state(1), state(2), direction) = (1 - lr) * qmatrix(state(1), state(2), direction) + lr *(rewardM + gamma * max(qmatrix(state1(1), state1(2), :)));

        for l = 1 : length(selectedRuns)
            if(i == selectedRuns(l))
                runSaver(:, l, cnt) = state;
            end
        end
        state = state1;
    end
    disp(i)
end

reward = [2, 12];
cat = [8, 4];
for ind = 1:length(selectedRuns)
    runStemp = zeros(2, 1000);
    for i = 1 : size(runStemp, 2)
        runStemp(:, i) = runSaver(:, ind, i);
    end

    runStemp(:, ~any(runStemp, 1)) = [];


    writerObj = VideoWriter(sprintf('%d.avi', selectedRuns(ind)));
    writerObj.FrameRate = 10;
    open(writerObj);
    figure
    for i = 1 : size(runStemp,2)
        map = 255 * ones(15, 15, 3);
        map(reward(1), reward(2), :) = [0 255 0];
        map(cat(1), cat(2), :) = [255 0 0];
        map(runStemp(1,i), runStemp(2,i), :) = [0 0 0];
        image(map);
        ax = gca;
        ax.YDir = 'normal';
        grid on
        title(sprintf("Trial %d", selectedRuns(ind)));
        pause(0.05)
        hold on
        writeVideo(writerObj, getframe(gcf));
        drawnow
    end
    close(writerObj);

    figure;
    plot(runStemp(2,:), runStemp(1,:), 'k', 'LineWidth', 1.5)
    ax = gca;
    ax.YDir = 'normal';
    xlim([1, 15])
    ylim([1, 15])
    hold on
    scatter(reward(2), reward(1), 'g', 'filled')
    scatter(cat(2), cat(1), 'r', 'filled')
    scatter(runStemp(2,1), runStemp(1,1), 'k', 'filled')
    title(sprintf("Path - Trial %d", selectedRuns(ind)))
    legend('Path - Rat', 'Reward', 'Cat', 'Start point')
    saveas(gcf, sprintf("Path_Trial%d.png", selectedRuns(ind)))
end

%% Question 2
figure
V_Y = qmatrix(:,:, 1) - qmatrix(:,:, 3);
V_X = -qmatrix(:,:, 2) + qmatrix(:,:, 4);
quiver(V_X, V_Y, 'b')
hold on
[M,c] = contour(sqrt(V_X.^2 + V_Y.^2));
c.LineWidth = 1.5;
title("Contour and gradients")
saveas(gcf, 'Contour_gradients.png')

%% Question 3
map = zeros(15, 15);
reward = [2, 12];
cat = [8, 4];
values = -1 * ones(15, 15);
values(sub2ind(size(values), [2, 8], [12, 4])) = [10, -10];
runNum = 10000;
eps = 0.1;
heatMap = zeros(21,20);

hm1 = 0;
for lr = 0.2:0.05:1.1
    hm1 = hm1 + 1;
    hm2 = 0;
    for gamma = 0.5:0.05:1.05
        hm2  = hm2 + 1;
        qmatrix = zeros(15, 15, 4);
        numberOfActions = zeros(1, runNum);
        i = 0;
        while i <= runNum
            i = i + 1;
            state = randi([1, 15], 1, 2);
            cnt = 0;
            while ~(isequal(state, reward) || isequal(state, cat))
                cnt = cnt + 1;
                if rand() < eps
                    direction = randi(4);
                else
                    [~, direction] = max(qmatrix(state(1), state(2), :));
                end

                if direction == 1
                    state1 = [state(1), min(state(2) + 1, 15)];
                elseif direction == 2
                    state1 = [state(1), max(state(2) - 1, 1)];
                elseif direction == 3
                    state1 = [max(state(1) - 1, 1), state(2)];
                elseif direction == 4
                    state1 = [min(state(1) + 1, 15), state(2)];
                end

                rewardM = values(state1(1), state1(2));
                qmatrix(state(1), state(2), direction) = (1 - lr) * qmatrix(state(1), state(2), direction) + lr *(rewardM + gamma * max(qmatrix(state1(1), state1(2), :)));

                state = state1;
            end
            numberOfActions(i) = cnt;
        end
        for s = 40:runNum
            if(mean(numberOfActions(s-39:s)) < 14)
                heatMap(hm1, hm2) = s;
                break;
            end
        end
    end
end

figure
heatmap(0.5:0.05:1.05,0.2:0.05:1.1, heatMap(1:19,1:12), 'Colormap', jet)
title("Heatmap")
xlabel("Discount factor")
ylabel("lr")
saveas(gcf, 'Heatmap.png')

%% Question 4
map = zeros(15, 15);
reward = [2, 12];
reward2 = [11, 11];
cat = [8, 4];
values = -1 * ones(15, 15);
values(sub2ind(size(values), [2, 8, 11], [12, 4, 11])) = [10, -10, 30];
qmatrix = zeros(15, 15, 4);
runNum = 10000;
eps = 0.1;
lr = 0.5;
gamma = 0.8;
runSaver = zeros(2, 11, 1000);
selectedRuns = [5, 212, 666];
i = 0;
while i <= runNum
    i = i + 1;
    state = randi([1, 15], 1, 2);
    cnt = 0;
    while ~(isequal(state, reward) || isequal(state, cat) || isequal(state, reward2))
        cnt = cnt + 1;
        if rand() < eps
            direction = randi(4);
        else
            [~, direction] = max(qmatrix(state(1), state(2), :));
        end

        if direction == 1
            state1 = [state(1), min(state(2) + 1, 15)];
        elseif direction == 2
            state1 = [state(1), max(state(2) - 1, 1)];
        elseif direction == 3
            state1 = [max(state(1) - 1, 1), state(2)];
        elseif direction == 4
            state1 = [min(state(1) + 1, 15), state(2)];
        end

        rewardM = values(state1(1), state1(2));
        qmatrix(state(1), state(2), direction) = (1 - lr) * qmatrix(state(1), state(2), direction) + lr *(rewardM + gamma * max(qmatrix(state1(1), state1(2), :)));

        for l = 1 : length(selectedRuns)
            if(i == selectedRuns(l))
                runSaver(:, l, cnt) = state;
            end
        end
        state = state1;
    end
    disp(i)
end

reward = [2, 12];
reward2 = [11, 11];
cat = [8, 4];
for ind = 1:length(selectedRuns)
    runStemp = zeros(2, 1000);
    for i = 1 : size(runStemp, 2)
        runStemp(:, i) = runSaver(:, ind, i);
    end

    runStemp(:, ~any(runStemp, 1)) = [];


    writerObj = VideoWriter(sprintf('Q4_%d.avi', selectedRuns(ind)));
    writerObj.FrameRate = 10;
    open(writerObj);
    figure
    for i = 1 : size(runStemp,2)
        map = 255 * ones(15, 15, 3);
        map(reward(1), reward(2), :) = [0 127 0];
        map(cat(1), cat(2), :) = [255 0 0];
        map(reward2(1), reward2(2), :) = [0 255 0];
        map(runStemp(1,i), runStemp(2,i), :) = [0 0 0];
        image(map);
        ax = gca;
        ax.YDir = 'normal';
        grid on
        title(sprintf("Trial %d", selectedRuns(ind)));
        pause(0.05)
        hold on
        writeVideo(writerObj, getframe(gcf));
        drawnow
    end
    close(writerObj);

    figure;
    plot(runStemp(2,:), runStemp(1,:), 'k', 'LineWidth', 1.5)
    ax = gca;
    ax.YDir = 'normal';
    xlim([1, 15])
    ylim([1, 15])
    hold on
    scatter(reward(2), reward(1), 'g', 'filled')
    alpha(.5)
    scatter(reward2(2), reward2(1), 'g', 'filled')
    scatter(cat(2), cat(1), 'r', 'filled')
    scatter(runStemp(2,1), runStemp(1,1), 'k', 'filled')
    title(sprintf("Path - Trial %d", selectedRuns(ind)))
    legend('Path - Rat', 'Reward', 'Reward2','Cat', 'Start point')
    saveas(gcf, sprintf("Q4_Path_Trial%d.png", selectedRuns(ind)))
end

figure
V_Y = qmatrix(:,:, 1) - qmatrix(:,:, 3);
V_X = -qmatrix(:,:, 2) + qmatrix(:,:, 4);
quiver(V_X, V_Y, 'b')
hold on
[M,c] = contour(sqrt(V_X.^2 + V_Y.^2));
c.LineWidth = 1.5;
title("Contour and gradients")
saveas(gcf, 'Q4_Contour_gradients.png')
%%
% Heatmap
map = zeros(15, 15);
reward = [2, 12];
reward2 = [11, 11];
cat = [8, 4];
values = -1 * ones(15, 15);
values(sub2ind(size(values), [2, 8, 11], [12, 4, 11])) = [10, -10, 30];
runNum = 10000;
eps = 0.1;
heatMap = zeros(21,20);

hm1 = 0;
for lr = 0.2:0.05:1.1
    hm1 = hm1 + 1;
    hm2 = 0;
    for gamma = 0.45:0.05:0.9
        hm2  = hm2 + 1;
        qmatrix = zeros(15, 15, 4);
        numberOfActions = zeros(1, runNum);
        i = 0;
        while i <= runNum
            i = i + 1;
            state = randi([1, 15], 1, 2);
            cnt = 0;
            while ~(isequal(state, reward) || isequal(state, cat) || isequal(state, reward2))
                cnt = cnt + 1;
                if rand() < eps
                    direction = randi(4);
                else
                    [~, direction] = max(qmatrix(state(1), state(2), :));
                end

                if direction == 1
                    state1 = [state(1), min(state(2) + 1, 15)];
                elseif direction == 2
                    state1 = [state(1), max(state(2) - 1, 1)];
                elseif direction == 3
                    state1 = [max(state(1) - 1, 1), state(2)];
                elseif direction == 4
                    state1 = [min(state(1) + 1, 15), state(2)];
                end

                rewardM = values(state1(1), state1(2));
                qmatrix(state(1), state(2), direction) = (1 - lr) * qmatrix(state(1), state(2), direction) + lr *(rewardM + gamma * max(qmatrix(state1(1), state1(2), :)));

                state = state1;
            end
            numberOfActions(i) = cnt;
        end
        for s = 40:runNum
            if(mean(numberOfActions(s-39:s)) < 14)
                heatMap(hm1, hm2) = s;
                break;
            end
        end
    end
end

figure
heatmap(0.45:0.05:0.9,0.2:0.05:1.1, heatMap(1:19,1:12), 'Colormap', jet)
title("Heatmap")
xlabel("Discount factor")
ylabel("lr")
saveas(gcf, 'Q4_Heatmap.png')

%% Question 5
map = zeros(15, 15);
reward = [2, 12];
reward2 = [11, 11];
cat = [8, 4];
values = -1 * ones(15, 15);
values(sub2ind(size(values), [2, 8, 11], [12, 4, 11])) = [10, -10, 30];
qmatrix = zeros(15, 15, 4);
runNum = 10000;
eps = 0.1;
lr = 0.5;
gamma = 0.8;
lambda = 0.8;
runSaver = zeros(2, 11, 1000);
selectedRuns = [5, 212, 666];
i = 0;
while i <= runNum
    i = i + 1;
    e = zeros(15, 15);
    state = randi([1, 15], 1, 2);
    cnt = 0;
    while ~(isequal(state, reward) || isequal(state, cat) || isequal(state, reward2))
        cnt = cnt + 1;
        if rand() < eps
            direction = randi(4);
        else
            [~, direction] = max(qmatrix(state(1), state(2), :));
        end

        if direction == 1
            state1 = [state(1), min(state(2) + 1, 15)];
        elseif direction == 2
            state1 = [state(1), max(state(2) - 1, 1)];
        elseif direction == 3
            state1 = [max(state(1) - 1, 1), state(2)];
        elseif direction == 4
            state1 = [min(state(1) + 1, 15), state(2)];
        end

        rewardM = values(state1(1), state1(2));
        e(state(1), state(2)) = e(state(1), state(2))+1;
        qmatrix(state(1), state(2), direction) = (1 - e(state(1), state(2)) * lr) * qmatrix(state(1), state(2), direction) + e(state(1), state(2)) * lr *(rewardM + gamma * max(qmatrix(state1(1), state1(2), :)));

        for l = 1 : length(selectedRuns)
            if(i == selectedRuns(l))
                runSaver(:, l, cnt) = state;
            end
        end
        e(state(1), state(2)) = e(state(1), state(2))* gamma * lambda;
        state = state1;
    end
    disp(i)
end

reward = [2, 12];
reward2 = [11, 11];
cat = [8, 4];
for ind = 1:length(selectedRuns)
    runStemp = zeros(2, 1000);
    for i = 1 : size(runStemp, 2)
        runStemp(:, i) = runSaver(:, ind, i);
    end

    runStemp(:, ~any(runStemp, 1)) = [];

    writerObj = VideoWriter(sprintf('Q5_%d.avi', selectedRuns(ind)));
    writerObj.FrameRate = 10;
    open(writerObj);
    figure
    for i = 1 : size(runStemp,2)
        map = 255 * ones(15, 15, 3);
        map(reward(1), reward(2), :) = [0 127 0];
        map(cat(1), cat(2), :) = [255 0 0];
        map(reward2(1), reward2(2), :) = [0 255 0];
        map(runStemp(1,i), runStemp(2,i), :) = [0 0 0];
        image(map);
        ax = gca;
        ax.YDir = 'normal';
        grid on
        title(sprintf("Trial %d", selectedRuns(ind)));
        pause(0.05)
        hold on
        writeVideo(writerObj, getframe(gcf));
        drawnow
    end
    close(writerObj);

    figure;
    plot(runStemp(2,:), runStemp(1,:), 'k', 'LineWidth', 1.5)
    ax = gca;
    ax.YDir = 'normal';
    xlim([1, 15])
    ylim([1, 15])
    hold on
    scatter(reward(2), reward(1), 'g', 'filled')
    alpha(.5)
    scatter(reward2(2), reward2(1), 'g', 'filled')
    scatter(cat(2), cat(1), 'r', 'filled')
    scatter(runStemp(2,1), runStemp(1,1), 'k', 'filled')
    title(sprintf("Path - Trial %d", selectedRuns(ind)))
    legend('Path - Rat', 'Reward', 'Reward2','Cat', 'Start point')
    saveas(gcf, sprintf("Q5_Path_Trial%d.png", selectedRuns(ind)))
end

figure
V_Y = qmatrix(:,:, 1) - qmatrix(:,:, 3);
V_X = -qmatrix(:,:, 2) + qmatrix(:,:, 4);
quiver(V_X, V_Y, 'b')
hold on
[M,c] = contour(sqrt(V_X.^2 + V_Y.^2));
c.LineWidth = 1.5;
title("Contour and gradients")
saveas(gcf, 'Q5_Contour_gradients.png')

% Heatmap
map = zeros(15, 15);
reward = [2, 12];
reward2 = [11, 11];
cat = [8, 4];
values = -1 * ones(15, 15);
values(sub2ind(size(values), [2, 8, 11], [12, 4, 11])) = [10, -10, 30];
runNum = 10000;
eps = 0.1;
heatMap = zeros(21,20);
lambda = 0.8;

hm1 = 0;
for lr = 0.2:0.05:1.1
    hm1 = hm1 + 1;
    hm2 = 0;
    for gamma = 0.45:0.05:0.9
        hm2  = hm2 + 1;
        qmatrix = zeros(15, 15, 4);
        numberOfActions = zeros(1, runNum);
        i = 0;
        while i <= runNum
            i = i + 1;
            e = zeros(15, 15);
            state = randi([1, 15], 1, 2);
            cnt = 0;
            while ~(isequal(state, reward) || isequal(state, cat) || isequal(state, reward2))
                cnt = cnt + 1;
                if rand() < eps
                    direction = randi(4);
                else
                    [~, direction] = max(qmatrix(state(1), state(2), :));
                end

                if direction == 1
                    state1 = [state(1), min(state(2) + 1, 15)];
                elseif direction == 2
                    state1 = [state(1), max(state(2) - 1, 1)];
                elseif direction == 3
                    state1 = [max(state(1) - 1, 1), state(2)];
                elseif direction == 4
                    state1 = [min(state(1) + 1, 15), state(2)];
                end

                rewardM = values(state1(1), state1(2));
                e(state(1), state(2)) = e(state(1), state(2))+1;
                qmatrix(state(1), state(2), direction) = (1 - e(state(1), state(2)) * lr) * qmatrix(state(1), state(2), direction) + e(state(1), state(2)) * lr *(rewardM + gamma * max(qmatrix(state1(1), state1(2), :)));
                e(state(1), state(2)) = e(state(1), state(2))* gamma * lambda;
                state = state1;
            end
            numberOfActions(i) = cnt;
        end
        for s = 40:runNum
            if(mean(numberOfActions(s-39:s)) < 14)
                heatMap(hm1, hm2) = s;
                break;
            end
        end
    end
end

figure
heatmap(0.45:0.05:0.9, 0.2:0.05:1.1, heatMap(1:19, 1:10), 'Colormap', jet)
title("Heatmap")
xlabel("Discount factor")
ylabel("lr")
saveas(gcf, 'Q5_Heatmap.png')

