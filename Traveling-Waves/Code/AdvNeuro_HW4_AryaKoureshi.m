clc; clear; close all;
%% Load data
load("data/ArrayData.mat")
load("data/CleanTrials.mat")
%% Initializing
fs = 200;
numChannels = numel(chan);

%% Removing bad trials
for ch = 1:numChannels
    chan(ch).lfp = chan(ch).lfp(:, Intersect_Clean_Trials);
end
numTrials = size(chan(1).lfp, 2);

%% LFP analysis
% Part a - Finding the most dominant frequency oscillation 
Ps = 0;
for ch = 1:numChannels
    LFPdata = zscore(chan(ch).lfp);
    for tr = 1:numTrials
        trial = LFPdata(:, tr);
        m = length(trial);
        n = pow2(nextpow2(m));
        Y = fftshift(fft(trial, n));
        Ps = Ps+abs(Y);
    end
end

normalizeConstant = log10((numTrials * numChannels) ^ 2 / n) * 10;

f = (-n/2:n/2-1) * (fs/n);
Ps_plot = log10(Ps.^2 / n) * 10;
pinkNoise = 1./f(n/2+2:end);
[~,~,spectrumRegressed] = regress(Ps_plot(n/2+2:end), pinkNoise');
pinkSpectrum = Ps_plot(n/2+2:end) - spectrumRegressed;

figure
loglog(f(n/2+2:end), pinkSpectrum, '--r')
hold on
loglog(f(n/2+2:end), Ps_plot(n/2+2:end), 'k')
title('Power spectrum obtained by averaging all trials across all channels (unnormalized)')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
xlim([0, 100])
legend('Estimated Pink Noise')

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Analysis_Part_a_1", '-dpng', '-r0')

figure
plot(f(n/2+2:end), pinkSpectrum, '--r')
hold on
plot(f(n/2+2:end), Ps_plot(n/2+2:end) - normalizeConstant, 'k')
title('Power spectrum obtained by averaging all trials across all channels (normalized)')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
xlim([0, 100])
legend('Estimated Pink Noise')

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Analysis_Part_a_2", '-dpng', '-r0')

% Removing pink noise
figure
hold on
cleanSpectrum = Ps_plot(n/2+2:end) - pinkSpectrum;
plot(f(n/2+2:end), Ps_plot(n/2+2:end)- normalizeConstant, 'r')
plot(f(n/2+2:end), cleanSpectrum - normalizeConstant, 'k')
legend('Original', 'Denoised (No Pink Noise)')
title('Power spectrum obtained by averaging all trials across all channels (normalized)')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
xlim([0, 100])
ylim([-5, 50])

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Analysis_Part_a_3", '-dpng', '-r0')


% Part b - Clustering electrodes based on their dominant oscillation frequency
figure
hold on
dominantFreq = nan * ChannelPosition;
normalizeConstant = log10((numTrials)^2) * 10;
for ch = 1:numChannels
    LFPdata = chan(ch).lfp;
    Ps = 0;
    for tr = 1:numTrials
        trialData = zscore(LFPdata(:, tr));
        m = length(trialData);
        n = pow2(nextpow2(m));
        Y = fftshift(fft(trialData, n));
        Ps = Ps + abs(Y);
    end

    f = (-n/2:n/2-1) * (fs/n);
    Ps = log10(Ps.^2/n) * 10;
    Ps_plot = removePinkNoise(Ps, f, n, 1);
    plot(f(n/2+2:end), Ps_plot(n/2+2:end) - normalizeConstant)
    [row, ~] = find(Ps_plot(n/2+2:end) == max(Ps_plot(n/2+2:end)));
    f_tmp = f(n/2+2:end);
    dmFreq = f_tmp(row);
    [row, col] = find(ChannelPosition==ch);
    dominantFreq(row, col) = dmFreq;
end

title('Power spectrum obtained by averaging all trials across all channels (normalized)')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
xlim([0, 70])
ylim([-35, 10])

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Analysis_Part_b_1", '-dpng', '-r0')

figure
plt = imagesc(dominantFreq);
set(plt,'AlphaData', ~isnan(dominantFreq))
colormap jet
c_bar = colorbar;
caxis([0, 14])
title('Dominant Frequencies')
ylabel(c_bar,'Frequency (Hz)')

textStrings = num2str(dominantFreq(:), '%0.2f');
textStrings = strtrim(cellstr(textStrings));  
[x, y] = meshgrid(1:10, 1:5);  
hStrings = text(x(:), y(:), textStrings(:), 'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  
textColors = repmat(dominantFreq(:) > midValue, 1, 3);  
set(hStrings, {'Color'}, num2cell(textColors, 2));

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Analysis_Part_b_2", '-dpng', '-r0')

% Part c - time-frequncy analysis of the LFP data
% STFT Method
STFTmap = 0;
for ch = 1:numChannels
    LFPdata = zscore(chan(ch).lfp);
    for tr = 1:numTrials
        trial = LFPdata(:, tr);
        [s,f,time_stft] = stft(trial, fs, 'Window', kaiser(60,5), 'OverlapLength', 40, 'FFTLength', fs);
        STFTmap = STFTmap + abs(s); 
    end
end

STFTmap_tmp = [];
for t = 1:size(STFTmap, 2)
    Ps = STFTmap(:, t);
    n = length(Ps);
    Ps_plot = removePinkNoise(Ps, f', n, 1);
    STFTmap_tmp(:, t) = Ps_plot;
end

figure
imagesc(time_stft - 1.2, f, STFTmap / (numTrials * numChannels));
title("STFT")
ylim([0, 40])
caxis([0, 14])
set(gca, 'YDir', 'normal')
c_bar = colorbar;
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
ylabel(c_bar, 'Power (dB)')

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Analysis_Part_c_1", '-dpng', '-r0')

figure
imagesc(time_stft - 1.2, f(n/2+2:end), STFTmap_tmp(n/2+2:end, :) / (numTrials * numChannels));
title("STFT + removed pink noise")
ylim([0, 40])
caxis([0, 14])
set(gca, 'YDir', 'normal')
c_bar = colorbar;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ylabel(c_bar, 'Power (dB)')

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Analysis_Part_c_2", '-dpng', '-r0')

% Welch Method
pxx_mean = 0;
for ch = 1:numChannels
    LFPdata = zscore(chan(ch).lfp);
    for tr = 1:numTrials
        trial = LFPdata(:, tr);
        tmp = buffer(trial, 80, 60);
        [pxx, f] = pwelch(tmp, 40, 20, fs, fs);
        pxx_mean = pxx_mean + pxx;
    end
end

pxx_clean = [];
for t = 1:size(pxx_mean, 2)
    Ps = pxx_mean(:, t);
    pxx_clean(:, t) = removePinkNoise(Ps, f', length(Ps), 2);
end

figure
imagesc(linspace(-1.2, 2, size(tmp, 2)), f, pxx_mean / (numTrials * numChannels))
c_bar = colorbar;
title("Power Spectrum (Welch)")
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
ylabel(c_bar, 'Power')
caxis([-0.05, 0.1])
set(gca, 'YDir', 'normal')
ylim([0, 40])

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Analysis_Part_c_3", '-dpng', '-r0')

figure
imagesc(linspace(-1.2, 2, size(tmp, 2)), f, pxx_clean / (numTrials * numChannels))
c_bar = colorbar;
title("Power Spectrum (Welch) + removed pink noise")
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
ylabel(c_bar, 'Power')
caxis([-0.05, 0.1])
set(gca, 'YDir', 'normal')
ylim([0, 40])

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Analysis_Part_c_4", '-dpng', '-r0')

% Wavelet
wt = 0;
for ch = 1:numChannels
    LFPdata = zscore(chan(ch).lfp);
    for tr = 1:numTrials
        trial = LFPdata(:, tr);
        [wt_tmp, f, ~] = cwt(trial,fs);
        wt = wt + abs(wt_tmp);
    end
end

figure
surface(Time, f, wt / (numTrials * numChannels));
axis tight
shading flat
ylim([1.2, 40])
title("Power Spectrum (Wavelet)")
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
c_bar = colorbar;
set(gca, 'YDir', 'normal')
ylabel(c_bar, 'Power (dB)')

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Analysis_Part_c_5", '-dpng', '-r0')

%% Phase propagation (Traveling waves)
% Part a - Bandpass filtering
dominantFreq = 12.5;
[b, a] = butter(3, [dominantFreq-1 dominantFreq+1] / (fs/2), 'bandpass');
freqz(b, a, fs, fs)
ax = findall(gcf, 'Type', 'axes');
set(ax, 'XLim', [10 20]);
xline(ax(1), dominantFreq, '--k');
xline(ax(2), dominantFreq, '--k');
title("The filter's frequency response")

for ch = 1:numChannels
    LFPdata = zscore(chan(ch).lfp);
    chan(ch).filtered_lfp = transpose(filtfilt(b, a, transpose(LFPdata)));
end

set(gcf, 'PaperPositionMode', 'auto')
print("report/Traveling_Waves_Part_a", '-dpng', '-r0')

% Part b - Calculating instantaneous pahse of filtered signals using HT
for ch = 1:numChannels
    instantaneous_phase = angle(hilbert(chan(ch).filtered_lfp));
    instantaneous_cos_phase = cos(angle(hilbert(chan(ch).filtered_lfp)));
    chan(ch).phase = instantaneous_phase;
    chan(ch).cos_phase = instantaneous_cos_phase;
end

% Part c - Displaying the traveling wave
clc; close all;
timesPlot = 1/fs:1/fs:3.2+1/fs;
timesPlot = floor(timesPlot*fs);
trial_no = 259;

finalAngle = zeros(size(ChannelPosition, 1), size(ChannelPosition, 2), length(timesPlot))*nan;
finalAngleCos = zeros(size(ChannelPosition, 1), size(ChannelPosition, 2), length(timesPlot))*nan;
finalData = zeros(size(ChannelPosition, 1), size(ChannelPosition, 2), length(timesPlot))*nan;
cntTime = 1;
for t = timesPlot
    for i = 1:size(ChannelPosition, 1)
        for j = 1:size(ChannelPosition, 2)
            ch_no = ChannelPosition(i, j);
            if ~isnan(ch_no)
                ang = chan(ch_no).cos_phase(:, trial_no);
                finalAngle(i, j, cntTime) = ang(cntTime);
                
                angCos = chan(ch_no).cos_phase(:, trial_no);
                finalAngleCos(i, j, cntTime) = angCos(cntTime);
                
                data = chan(ch_no).filtered_lfp(:, trial_no);
                data = data(cntTime);
                finalData(i, j, cntTime) = ang(cntTime);
            end
        end
    end
    cntTime = cntTime+1;
end

focus_chs = [11, 16, 21, 26, 31, 36];

fig = figure('units','normalized','outerposition',[0 0 1 1]);
writerObj = VideoWriter("report/videos/Trial_" + num2str(trial_no));
writerObj.FrameRate = 25;
open(writerObj);

cntTime = 101;
phi_0 = 0;
for t = timesPlot(cntTime:end-cntTime)
    time_2_plot = t/fs-1/fs-1.2;
    subplot(2,2,[1, 2])
    indexes = cntTime-100:cntTime+100;
    times_2_plot = indexes/fs-1/fs-1.2;
    for ch_no = 1:48
        [row, col] = find(ChannelPosition==ch_no);
        frame_data_2_plot = finalAngle(row, col, indexes);
        frame_data_2_plot = reshape(frame_data_2_plot, [], length(indexes));
        plot_color = [0.9 0.9 0.9];
        plot_line_width = 1;
        plot(times_2_plot, frame_data_2_plot, 'color', plot_color, 'LineWidth', plot_line_width)
        hold on
    end
    
    ch_counter = 1;
    plot_color = [0, 0, 0];
    for ch_no = focus_chs
        [row, col] = find(ChannelPosition==ch_no);
        frame_data_2_plot = finalAngle(row, col, indexes);
        frame_data_2_plot = reshape(frame_data_2_plot, [], length(indexes));
        plot_color = plot_color + 0.1;
        plot_line_width = 2;
        ch_counter = ch_counter+1;
        plot(times_2_plot, frame_data_2_plot, 'color', plot_color, 'LineWidth', plot_line_width)
    end
    title("Trial No."+num2str(trial_no))
    xline(time_2_plot, '--r');
    hold off
    ylim([-1.5 , 1.5])
    xlim([times_2_plot(1), times_2_plot(end)])
    xlabel('Time (s)')
    
    
    subplot(2,2,3)
    frame_angle_2_plot = finalAngle(:, :, cntTime);
    img = imagesc(frame_angle_2_plot);
    colormap hot
    set(img,'AlphaData', ~isnan(frame_angle_2_plot))
    hold on
    ch_counter = 1;
    plot_color = [0, 0, 0];
    for ch_no = focus_chs
        [row, col] = find(ChannelPosition==ch_no);
        plot_color = plot_color + 0.1;
        ch_counter = ch_counter+1;
        scatter(col, row, 80, plot_color, 'filled')
    end
    title("Time = "+num2str(time_2_plot))
    hold off    
    colorbar
    caxis([-1, 1])
    
    subplot(2,2,4)
    [u,v] = gradient(frame_angle_2_plot,1,1);
    [y,x] = ndgrid(1:5,1:10);
    quiver(x,y,u,v)
    hold on
    
    u_mean = nanmean(u, 'all');
    v_mean = nanmean(v, 'all');
    u_mean = 10*u_mean;
    v_mean = 10*v_mean;
    
    u = u/400e-6;
    v = v/400e-6;
    
    plot([5, 5+u_mean], [3, 3+v_mean], 'k->', 'LineWidth', 2)
    
    ch_counter = 1;
    plot_color = [0, 0, 0];
    for ch_no = focus_chs
        [row, col] = find(ChannelPosition==ch_no);
        plot_color = plot_color + 0.1;
        ch_counter = ch_counter+1;
        scatter(col, row, 80, plot_color, 'filled')
    end
    
    ylim([0.5, 5.5])
    xlim([0.5, 10.5])
    hold off
    pgd = pgdCalculator(u, v);
    d_phi = frame_angle_2_plot - phi_0;
    speed = 100*speedCalculator(u, v, d_phi, 1/fs);
    phi_0 = frame_angle_2_plot;
    title("PGD = "+num2str(pgd)+" | Speed = "+num2str(speed))
    legend('Gradient Vectors', 'Avg Gradient Vector')
    
    axis ij

    cntTime = cntTime+1;
    frame = getframe(fig);
    for frame_index = 1:4
        writeVideo(writerObj,frame);
    end
end
close(writerObj)


% Calculating parameters for all trials
timesPlot = 1/fs:1/fs:3.2+1/fs;
timesPlot = floor(timesPlot*fs);
numTrials = size(chan(1).lfp, 2);
finalAngle = zeros(size(ChannelPosition, 1), size(ChannelPosition, 2), length(timesPlot)) * nan;
finalPDG = zeros(numTrials, length(timesPlot));
finalSpeed = zeros(numTrials, length(timesPlot));
finalGradientDirectionMean = zeros(numTrials, length(timesPlot));
finalGradientDirectionAll = zeros(numTrials, length(timesPlot), 5, 5)*nan;

for tr = 1:size(chan(1).lfp, 2)
    cntTime = 1;
    phi_0 = 0;
    for t = timesPlot
        for i = 1:size(ChannelPosition, 1)
            for j = 1:size(ChannelPosition, 2)
                ch = ChannelPosition(i, j);
                if ~isnan(ch)
                    ang = chan(ch).phase(:, tr);
                    finalAngle(i, j, cntTime) = ang(cntTime);
                end
            end
        end
        [u, v] = gradient(finalAngle(:,:,cntTime),1,1);
        u = u/400e-6;
        v = v/400e-6;
        finalPDG(tr, cntTime) = pgdCalculator(u, v);
        d_phi = finalAngle(:,:,cntTime) - phi_0;
        finalSpeed(tr, cntTime) = speedCalculator(u, v, d_phi, 1/fs);
        phi_0 = finalAngle(:,:,cntTime);
        for row_grad = 1:size(u, 1)
            for col_grad = 1:size(u, 2)
                ch = ChannelPosition(row_grad, col_grad);
                if ~isnan(ch)
                    a = [u(row_grad, col_grad), v(row_grad, col_grad), 0];
                    b = [1, 0, 0];
                    finalGradientDirectionAll(tr, cntTime, row_grad, col_grad) = atan2(norm(cross(a,b)),dot(a,b))*180/pi;
                end
            end
        end
        a = [nanmean(u, 'all'), nanmean(v, 'all'), 0];
        b = [1, 0, 0];
        finalGradientDirectionMean(tr, cntTime) = atan2(norm(cross(a,b)),dot(a,b))*180/pi;
        cntTime = cntTime+1;
    end
end

% Plotting hitogram of PGDs and gradient directions
figure
histogram(finalGradientDirectionMean(:,1:241), 200,'Normalization', 'pdf')
hold on
histogram(finalGradientDirectionMean(:,241:end), 200,'Normalization', 'pdf')
title('PDF of Direction of Wave Propagation')
legend('Before Onset', 'After Onset')
xlabel('Propagation Direction (degree)')
ylabel('Probability Density')

set(gcf, 'PaperPositionMode', 'auto')
print("report/Traveling_Waves_Part_e_1", '-dpng', '-r0')

figure
finalGradientDirectionAll_tmp = finalGradientDirectionAll(:, 1:241, :, :);
finalGradientDirectionAll_tmp = finalGradientDirectionAll_tmp(~isnan(finalGradientDirectionAll_tmp));
histogram(finalGradientDirectionAll_tmp, 200,'Normalization', 'pdf')
hold on
finalGradientDirectionAll_tmp = finalGradientDirectionAll(:, 241:end, : ,:);
finalGradientDirectionAll_tmp = finalGradientDirectionAll_tmp(~isnan(finalGradientDirectionAll_tmp));
histogram(finalGradientDirectionAll_tmp, 200,'Normalization', 'pdf')
title('PDF of Direction of Gradient of all channels')
legend('Before Onset', 'After Onset')
xlabel('Gradient Direction (degree)')
ylabel('Probability Density')

set(gcf, 'PaperPositionMode', 'auto')
print("report/Traveling_Waves_Part_e_2", '-dpng', '-r0')

figure
histogram(finalSpeed(:,1:241)*100, 200,'Normalization', 'pdf')
hold on
histogram(finalSpeed(:,241:end)*100, 200,'Normalization', 'pdf')
title('PDF of Speed of Wave Propagation')
legend('Before Onset', 'After Onset')
xlabel('Speed (cm/s)')
ylabel('Probability Density')
xlim([0, 50])

set(gcf, 'PaperPositionMode', 'auto')
print("report/Traveling_Waves_Part_e_3", '-dpng', '-r0')

% Repeating last part just for times with pdg > 0.5

[rows, cols] = find(finalPDG<=0.5);

finalGradientDirectionMeanBackup = finalGradientDirectionMean;
finalGradientDirectionAllBackup = finalGradientDirectionAll;
finalSpeedBackup = finalSpeed;

for i = 1:length(rows)
    row = rows(i);
    col = cols(i);
    finalGradientDirectionMeanBackup(row,col) = nan;
    finalGradientDirectionAllBackup(row,col, :, :) = nan;
    finalSpeedBackup(row,col) = nan;
end

figure
histogram(finalGradientDirectionMeanBackup, 200,'Normalization', 'pdf')
title('PDF of Direction of Wave Propagation')
xlabel('Propagation Direction (degree)')
ylabel('Probability Density')

set(gcf, 'PaperPositionMode', 'auto')
print("report/Traveling_Waves_Part_e_4", '-dpng', '-r0')

figure
finalGradientDirectionAll_tmp = finalGradientDirectionAllBackup;
finalGradientDirectionAll_tmp = finalGradientDirectionAll_tmp(~isnan(finalGradientDirectionAll_tmp));
histogram(finalGradientDirectionAll_tmp, 200,'Normalization', 'pdf')
title('PDF of Direction of Gradient of all channels')
xlabel('Gradient Direction (degree)')
ylabel('Probability Density')

set(gcf, 'PaperPositionMode', 'auto')
print("report/Traveling_Waves_Part_e_5", '-dpng', '-r0')

figure
histogram(finalSpeedBackup*100, 200, 'Normalization', 'pdf')
title('PDF of Speed of Wave Propagation')
xlabel('Speed (cm/s)')
ylabel('Probability Density')
xlim([0, 50])

set(gcf, 'PaperPositionMode', 'auto')
print("report/Traveling_Waves_Part_e_6", '-dpng', '-r0')

% Part f - Validating the observations
figure
h = histogram(finalGradientDirectionMeanBackup, 200, 'Normalization', 'pdf');
values = h.Values;
bins = h.BinEdges(1:end-1)-90;
values = [values(end/2+1:end) values(1:end/2)];
plot(bins,values)
xlabel('Propagation Direction (degree)')
ylabel('Probability Density')
title('PDF of Direction of Wave Propagation (wraped between -90 and 90)')

set(gcf, 'PaperPositionMode', 'auto')
print("report/Traveling_Waves_Part_f", '-dpng', '-r0')

%% Functions

function Ps_plot = removePinkNoise(Ps, f, n, bool)
    if bool == 1
        pinkNoise = 1./f(n/2+2:end);
        [~,~,spectrumRegressed] = regress(Ps(n/2+2:end), pinkNoise');
        pinkSpectrum = Ps(n/2+2:end) - spectrumRegressed;
        Ps_plot = Ps;
        Ps_plot(n/2+2:end) = Ps(n/2+2:end) - pinkSpectrum;
    else
        pinkNoise = 1./f;
        if pinkNoise(1)==Inf, pinkNoise(1)=pinkNoise(2); end
        [~,~,spectrumRegressed] = regress(Ps, pinkNoise');
        pinkSpectrum = Ps-spectrumRegressed;
        Ps_plot = Ps;
        Ps_plot = Ps-pinkSpectrum;
    end
end

function pgd = pgdCalculator(fx, fy)
    pgd = norm(nanmean(fx, 'all'), nanmean(fy, 'all')) / nanmean(sqrt(fx.^2+fy.^2), 'all');
end

function speed = speedCalculator(fx, fy, dphi, dt)
    speed = abs(nanmean(dphi/dt, 'all')) / nanmean(sqrt(fx.^2+fy.^2), 'all');
end