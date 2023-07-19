clc; clear; close all;

load("data/ArrayData.mat")
load("data/CleanTrials.mat")

fs = 200;
numChannels = numel(chan);
numTrials = size(Intersect_Clean_Trials, 1);

%% Part a
Ps = 0;
for ch = 1:numChannels
    LFPdata = zscore(chan(ch).lfp(:, Intersect_Clean_Trials));
    for tr = 1:numTrials
        trial = LFPdata(:, tr);
        n = pow2(nextpow2(length(trial)));
        Ps = Ps + abs(fftshift(fft(trial, n)));
    end
end

f = linspace(-fs/2, fs/2 - fs/n, n);
f = f((n/2+2):end);

Ps_plot = log10(Ps.^2 / n) * 10;
pinkNoise = 1./f;
[~,~,spectrumRegressed] = regress(Ps_plot(n/2+2:end), pinkNoise');
pinkSpectrum = Ps_plot(n/2+2:end) - spectrumRegressed;

figure
hold on
plot(f, pinkSpectrum, ':b')
cleanSpectrum = Ps_plot(n/2+2:end) - pinkSpectrum;
plot(f, Ps_plot(n/2+2:end)- log10((numTrials * numChannels) ^ 2 / n) * 10, 'r', LineWidth=1.5)
plot(f, cleanSpectrum - log10((numTrials * numChannels) ^ 2 / n) * 10, 'g', LineWidth=1.5)
legend('Pink Noise', 'Original', 'Denoised (No Pink Noise)')
title('Average over all trials and all channels (normalized)')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
ylim([-10, 60])

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Analysis_Part_a", '-dpng', '-r0')

%% Part b
dominantFreq = nan * ChannelPosition;
for ch = 1:numChannels
    LFPdata = chan(ch).lfp;
    Ps = 0;
    for tr = 1:numTrials
        trialData = zscore(LFPdata(:, tr));
        n = pow2(nextpow2(length(trialData)));
        Ps = Ps + abs(fftshift(fft(trialData, pow2(nextpow2(length(trialData))))));
    end
    f = linspace(-fs/2, fs/2 - fs/n, n);
    f = f((n/2+2):end);

    Ps = log10(Ps.^2/n) * 10;
    Ps = Ps(n/2+2:end);
    pinkNoise = 1./f;
    if pinkNoise(1)==Inf, pinkNoise(1)=pinkNoise(2); end
    [~,~,spectrumRegressed] = regress(Ps, pinkNoise');
    Ps_plot = spectrumRegressed;
    [row, ~] = find(Ps_plot == max(Ps_plot));
    [row1, col1] = find(ChannelPosition==ch);
    dominantFreq(row1, col1) = f(row);
end

figure
plt = imagesc(dominantFreq);
colormap parula
caxis([0, 13])
title('Dominant Frequencies - Electrodes')
ylabel(colorbar,'Frequency (Hz)')
textStrings = num2str(dominantFreq(:), '%0.2f');
textStrings = strtrim(cellstr(textStrings));  
[x, y] = meshgrid(1:10, 1:5);  
hStrings = text(x(:), y(:), textStrings(:), 'HorizontalAlignment', 'center');

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Analysis_Part_b", '-dpng', '-r0')

%% Part c
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
    pinkNoise = 1./f(n/2+2:end)';
    [~,~,spectrumRegressed] = regress(Ps(n/2+2:end), pinkNoise');
    pinkSpectrum = Ps(n/2+2:end) - spectrumRegressed;
    Ps_plot = Ps;
    Ps_plot(n/2+2:end) = spectrumRegressed;
    STFTmap_tmp(:, t) = Ps_plot;
end

figure
imagesc(time_stft - 1.2, f, STFTmap / (numTrials * numChannels));
title("STFT")
ylim([0, 100])
caxis([0, 14])
set(gca, 'YDir', 'normal')
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
ylabel(colorbar, 'Power (dB)')

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Analysis_Part_c_1", '-dpng', '-r0')

figure
imagesc(time_stft - 1.2, f, STFTmap_tmp / (numTrials * numChannels));
title("STFT (Removed pink noise)")
ylim([0, 100])
caxis([0, 14])
set(gca, 'YDir', 'normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ylabel(colorbar, 'Power (dB)')

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Analysis_Part_c_2", '-dpng', '-r0')

%% Phase propagation (Traveling waves)
%% Part a
[b, a] = butter(3, [11.5 13.5] / (fs/2), 'bandpass');
freqz(b, a, fs, fs)
title("Butterworth filter - order=3")
set(gcf, 'PaperPositionMode', 'auto')
print("report/Traveling_Waves_Part_a", '-dpng', '-r0')

%% Part b
Ps = 0;
for ch = 1:numChannels
    LFPdata = zscore(chan(ch).lfp);
    chan(ch).filtered_lfp = transpose(filtfilt(b, a, transpose(LFPdata)));
    LFPdata = zscore(chan(ch).filtered_lfp(:, Intersect_Clean_Trials));
    for tr = 1:numTrials
        trial = LFPdata(:, tr);
        n = pow2(nextpow2(length(trial)));
        Ps = Ps + abs(fftshift(fft(trial, n)));
    end
    instantaneous_phase = angle(hilbert(chan(ch).filtered_lfp));
    instantaneous_cos_phase = cos(angle(hilbert(chan(ch).filtered_lfp)));
    chan(ch).phase = instantaneous_phase;
    chan(ch).cos_phase = instantaneous_cos_phase;
end

f = linspace(-fs/2, fs/2 - fs/n, n);
f = f((n/2+2):end);

Ps_plot = log10(Ps.^2 / n) * 10;
pinkNoise = 1./f;
[~,~,spectrumRegressed] = regress(Ps_plot(n/2+2:end), pinkNoise');
pinkSpectrum = Ps_plot(n/2+2:end) - spectrumRegressed;

figure
hold on
plot(f, pinkSpectrum, ':b')
cleanSpectrum = Ps_plot(n/2+2:end) - pinkSpectrum;
plot(f, Ps_plot(n/2+2:end)- log10((numTrials * numChannels) ^ 2 / n) * 10, 'r', LineWidth=1.5)
plot(f, cleanSpectrum - log10((numTrials * numChannels) ^ 2 / n) * 10, 'g', LineWidth=1.5)
legend('Pink Noise', 'Filtered LFP', 'Denoised Filtered LFP (No Pink Noise)')
title('Average over all trials and all channels (normalized)')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
ylim([-10, 60])

set(gcf, 'PaperPositionMode', 'auto')
print("report/LFP_Traveling_Waves_Part_b", '-dpng', '-r0')

%% Part c
clc; close all;
timesPlot = floor(1:1:3.2*fs+1);
trial_no = 259;
fig = figure('units','normalized','outerposition',[0 0 1 1]);
writerObj = VideoWriter("report/videos/Trial_" + num2str(trial_no));
writerObj.FrameRate = 25;
open(writerObj);

finalAngleCos = zeros(size(ChannelPosition, 1), size(ChannelPosition, 2), length(timesPlot)) * nan;
cntTime = 1;
phi_0 = 0;
for t = timesPlot
    for i = 1:size(ChannelPosition, 1)
        for j = 1:size(ChannelPosition, 2)
            if ~isnan(ChannelPosition(i, j))       
                angCos = chan(ChannelPosition(i, j)).cos_phase(:, trial_no);
                finalAngleCos(i, j, cntTime) = angCos(cntTime);
            end
        end
    end

    subplot(2,1,1)
    frame_angle_2_plot = finalAngleCos(:, :, cntTime);
    img = imagesc(finalAngleCos(:, :, cntTime));
    colormap parula
    hold on

    title("Time: " + num2str((t-1)/fs))
    hold off    
    colorbar
    caxis([-1, 1])
    
    subplot(2,1,2)
    [u,v] = gradient(finalAngleCos(:, :, cntTime),1,1);
    [y,x] = ndgrid(1:5,1:10);
    quiver(x,y,u,v,'r')
    hold on
    
    u = u/400e-6;
    v = v/400e-6;
        
    ylim([0.5, 5.5])
    xlim([0.5, 10.5])
    hold off
    pgd = norm(nanmean(u, 'all'), nanmean(v, 'all')) / nanmean(sqrt(u.^2+v.^2), 'all');
    d_phi = finalAngleCos(:, :, cntTime) - phi_0;
    speed = 100 * abs(nanmean(d_phi * fs, 'all')) / nanmean(sqrt(u.^2+v.^2), 'all');
    phi_0 = finalAngleCos(:, :, cntTime);
    title("PGD = " + num2str(pgd) + ", Speed = " + num2str(speed))
    legend('Gradient Vectors')
    
    axis ij

    frame = getframe(fig);
    for frame_index = 1:2
        writeVideo(writerObj,frame);
    end

    cntTime = cntTime+1;
end
close(writerObj)

finalAngleCos = zeros(size(ChannelPosition, 1), size(ChannelPosition, 2), length(timesPlot)) * nan;
finalPDG = zeros(numTrials, length(timesPlot));
finalSpeed = zeros(numTrials, length(timesPlot));
finalGradientDirectionAll = zeros(numTrials, length(timesPlot), 5, 5) * nan;
finalGradientDirectionMean = zeros(numTrials, length(timesPlot), 5, 5) * nan;
for tr = 1:numTrials
    cntTime = 1;
    phi_0 = 0;
    for t = timesPlot
        for i = 1:size(ChannelPosition, 1)
            for j = 1:size(ChannelPosition, 2)
                if ~isnan(ChannelPosition(i, j))       
                    angCos = chan(ChannelPosition(i, j)).cos_phase(:, trial_no);
                    finalAngleCos(i, j, cntTime) = angCos(cntTime);
                end
            end
        end

        [u,v] = gradient(finalAngleCos(:, :, cntTime),1,1);      
        u = u/400e-6;
        v = v/400e-6;

        for i = 1:size(ChannelPosition, 1)
            for j = 1:size(ChannelPosition, 2)
                if ~isnan(ChannelPosition(i, j))
                    a = [u(i, j), v(i, j), 0];
                    b = [1, 0, 0];
                    finalGradientDirectionAll(tr, cntTime, i, j) = atan2(norm(cross([u(i, j), v(i, j), 0], [1, 0, 0])), dot([u(i, j), v(i, j), 0], [1, 0, 0])) * 180 / pi;
                end
            end
        end
        finalGradientDirectionMean(tr, cntTime) = atan2(norm(cross([nanmean(u, 'all'), nanmean(v, 'all'), 0], [1, 0, 0])),dot([nanmean(u, 'all'), nanmean(v, 'all'), 0], [1, 0, 0])) * 180 / pi;

        finalPDG(tr, cntTime) = norm(nanmean(u, 'all'), nanmean(v, 'all')) / nanmean(sqrt(u.^2+v.^2), 'all');
        d_phi = finalAngleCos(:, :, cntTime) - phi_0;
        finalSpeed(tr, cntTime) = abs(nanmean(d_phi * fs, 'all')) / nanmean(sqrt(u.^2+v.^2), 'all');
        phi_0 = finalAngleCos(:, :, cntTime);
        cntTime = cntTime+1;
    end
end

%% Part e
figure
finalGradientDirectionAll = finalGradientDirectionAll(~isnan(finalGradientDirectionAll(:, :, :, :)));
histogram(finalGradientDirectionAll, 400,'Normalization', 'pdf', 'FaceColor', 'g')
title('Direction of Gradients over all channels')
xlabel('Gradient Direction (degree)')
ylabel('Probability Density')
xlim([0, 180])
set(gcf, 'PaperPositionMode', 'auto')
print("report/Traveling_Waves_Part_e_1", '-dpng', '-r0')

figure
histogram(finalSpeed*100, 400,'Normalization', 'pdf', 'FaceColor', 'r')
title('Speed')
xlabel('Speed (cm/s)')
ylabel('Probability Density')
xlim([0, 50])

set(gcf, 'PaperPositionMode', 'auto')
print("report/Traveling_Waves_Part_e_2", '-dpng', '-r0')

%% Part f
[r, c] = find(finalPDG <= 0.5);
for i = 1:length(r)
    finalGradientDirectionMean(r(i),c(i)) = nan;
end

figure
histogram(finalGradientDirectionMean, 400, 'Normalization', 'pdf', 'FaceColor', 'k');
xlabel('Propagation Direction (degree)')
ylabel('Probability Density')
title('PDF of Direction of Wave Propagation')
set(gcf, 'PaperPositionMode', 'auto')
print("report/Traveling_Waves_Part_f", '-dpng', '-r0')
