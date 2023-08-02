clc; clear; close all;
%% Add to path
addpath('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/DatabaseCode/DatabaseCode/');
addpath('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/SaliencyToolbox');
addpath('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Our saliency model/JuddSaliencyModel/JuddSaliencyModel');
addpath('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/matlabPyrTools-master');
addpath('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/LabelMeToolbox-master/');
addpath('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/LabelMeToolbox-master/features');
addpath('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/LabelMeToolbox-master/imagemanipulation');
addpath('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/FaceDetect');
addpath('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/voc-release5');
addpath('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/gbvs/gbvs')
addpath('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/gbvs/gbvs/util');

%% Eye tracking database - One user
showEyeData('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/DATA/DATA/hp', ...
            'D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/ALLSTIMULI/ALLSTIMULI') % one user

%% Eye tracking database - All users
showEyeDataAcrossUsers('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/ALLSTIMULI/ALLSTIMULI', 3) % all users

%%  Saliency model
pics = [100, 200, 300, 400, 500, 600];
filePattern = fullfile('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/DATA/DATA/hp', '*.mat');
files = dir(filePattern);
fileNames = cell(1, size(files, 1));
[fileNames{:}] = deal(files.name);

for pic = 1:size(pics, 2)
    imagefile = fileNames{pics(pic)};
    features = saliency("D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/ALLSTIMULI/ALLSTIMULI/" + imagefile(1:end-4) + ".jpeg");
    
    load model;
    data = load(fullfile("D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/DATA/DATA/hp", imagefile));
    data = data.(imagefile(1:end-4)).DATA.eyeData;
    
    XFirstHalf = data(1:floor(end/2), 1);
    YFirstHalf = data(1:floor(end/2), 2);
    
    X_second_half = data(floor(end/2)+1:end, 1);
    Y_second_half = data(floor(end/2)+1:end, 2);
    
    img = imread("D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/ALLSTIMULI/ALLSTIMULI/" + imagefile(1:end-4) + ".jpeg");
    [w_img, h, c] = size(img);
    dims = [200, 200];
    
    figure('Visible', 'off');
    imshow(img);
    hold on;
    plot(XFirstHalf, YFirstHalf, '.y');
    scatter(XFirstHalf(1), YFirstHalf(1), 'b', 'filled');
    scatter(XFirstHalf(end), YFirstHalf(end), 'r', 'filled');
    
    set(gcf, 'PaperPositionMode', 'auto')
    print("D:/Sharif/Term2/Advanced Neuroscience/HW/8/images/saliency/part1_" + num2str(pics(pic)), '-dpng', '-r0')
    
    meanVec = model.whiteningParams(1, :);
    stdVec = model.whiteningParams(2, :);
    features = (features - meanVec) ./ stdVec;
    
    for type = 1:8
        w = model.w;
        
        switch type
            case 1 % Subband
                w(:, 1:12) = 0;
                w(:, 14:end) = 0;
            case 2 % Itti
                w(:, 1:13) = 0;
                w(:, 17:end) = 0;
            case 3 % Color
                w(1:16) = 0;
                w(:, 28:end) = 0;
            case 4 % Torralba
                w(:, 1:27) = 0;
                w(:, 29:end) = 0;
            case 5 % Horizon
                w(:, 1:28) = 0;
                w(:, 30:end) = 0;
            case 6 % Object
                w(:, 1:29) = 0;
                w(:, 32:end) = 0;
            case 7 % Center
                w(:, 1:32) = 0;
            case 8 % All
                % No modification needed
        end
    
        saliencyMap = (features * w(1:end-1)') + w(end);
        saliencyMap = (saliencyMap - min(saliencyMap)) / (max(saliencyMap) - min(saliencyMap));
        saliencyMap = reshape(saliencyMap, dims);
        saliencyMap = imresize(saliencyMap, [w_img, h]);
    
        figure('Visible', 'off');
        imshow(saliencyMap * 255 / max(saliencyMap, [], 'all'));
        colormap gray
        hold on;
        plot(XFirstHalf, YFirstHalf, '.y');
        scatter(XFirstHalf(1), YFirstHalf(1), 'blue', 'filled');
        scatter(XFirstHalf(end), YFirstHalf(end), 'red', 'filled');
    
        set(gcf, 'PaperPositionMode', 'auto')
        print("D:/Sharif/Term2/Advanced Neuroscience/HW/8/images/saliency/part1_" + num2str(pics(pic)) + "_type" + num2str(type), '-dpng', '-r0')
    end
end

%% Eye positoin on saliency map
load model;

filePattern = fullfile('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/DATA/DATA/hp', '*.mat');
files = dir(filePattern);
fileNames = cell(1, size(files, 1));
[fileNames{:}] = deal(files.name);

subjectNames = ["CNG", "ajs", "emb", "ems", "ff", "hp", "jcw", "jw", "kae", "krl", "po", "tmj", "tu", "ya", "zb"];

for subjectName = 1:size(subjectNames, 2)
    fileCounter = 1;
    for file = 1:size(fileNames, 2)
        fileName = fileNames(file);
        fileName = fileName{1}(1:end-4);
    
        dataFilePath = fullfile("D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/DATA/DATA", subjectNames(subjectName), fileName + ".mat");
        eyeData = load(dataFilePath);
        eyeData = eyeData.(fileName).DATA.eyeData;
    
        imgPath = fullfile("D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/ALLSTIMULI/ALLSTIMULI", fileName + ".jpeg");
        
        img = imread(imgPath);
        [imgWidth, imgHeight, ~] = size(img);
        resizedDims = [200, 200];
        map{1} = img;

        features = saliency(imgPath);
        meanVec = model.whiteningParams(1, :);
        stdVec = model.whiteningParams(2, :);
        features = (features - meanVec) ./ stdVec;

        saliencyMap = (features * model.w(1:end-1)') + model.w(end);
        saliencyMap = (saliencyMap - min(saliencyMap)) / (max(saliencyMap) - min(saliencyMap));
        saliencyMap = reshape(saliencyMap, resizedDims);
        saliencyMap = imresize(saliencyMap, [imgWidth, imgHeight]);

        XFirstHalf = eyeData(1:floor(end/2), 1);
        YFirstHalf = eyeData(1:floor(end/2), 2);
        figure;
        imshow(255 * saliencyMap / max(saliencyMap, [], 'all'));
        hold on;
        plot(XFirstHalf, YFirstHalf, '.y');
        scatter(XFirstHalf(1), YFirstHalf(1), 'b', 'filled');
        scatter(XFirstHalf(end), YFirstHalf(end), 'r', 'filled');
        colormap gray;

        set(gcf, 'PaperPositionMode', 'auto')
        print("D:/Sharif/Term2/Advanced Neuroscience/HW/8/images/saliency/part2_sub" + num2str(subjectNames(subjectName)) + "_" + num2str(fileName) + "_saliency", '-dpng', '-r0')  

        close all;
        if fileCounter == 5
            break;
        end
        fileCounter = fileCounter + 1;
    end
end

%%  You want to compare saliency maps to fixations
load model;

filePattern = fullfile('D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/DATA/DATA/hp', '*.mat');
files = dir(filePattern);
fileNames = cell(1, size(files, 1));
[fileNames{:}] = deal(files.name);

subjectNames = ["CNG", "ajs", "emb", "ems", "ff", "hp", "jcw", "jw", "kae", "krl", "po", "tmj", "tu", "ya", "zb"];
ROCs = zeros(size(subjectNames, 2), 2, 8, 2);

for subjectName = 1:size(subjectNames, 2)
    fileCounter = 1;
    for file = 1:size(fileNames, 2)
        fileName = fileNames(file);
        fileName = fileName{1}(1:end-4);
    
        dataFilePath = fullfile("D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/DATA/DATA", subjectNames(subjectName), fileName + ".mat");
        eyeData = load(dataFilePath);
        eyeData = eyeData.(fileName).DATA.eyeData;
    
        imgPath = fullfile("D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/ALLSTIMULI/ALLSTIMULI", fileName + ".jpeg");
        fixationMapPath = fullfile("D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Eye tracking database/ALLFIXATIONMAPS/ALLFIXATIONMAPS", fileName + "_fixMap.jpg");
    
        saliencyMapPath = fullfile("D:/Sharif/Term2/Advanced Neuroscience/HW/8/files/Our saliency model/SaliencyMaps/SaliencyMaps", fileName + "SM.jpg");
    
        img = imread(imgPath);
        [imgWidth, imgHeight, ~] = size(img);
        resizedDims = [200, 200];
        map{1} = img;
    
        features = saliency(imgPath);
        meanVec = model.whiteningParams(1, :);
        stdVec = model.whiteningParams(2, :);
        features = (features - meanVec) ./ stdVec;
        
        for type = 1:8
            w = model.w;
            
            switch type
                case 1 % Subband
                    w(:, 1:12) = 0;
                    w(:, 14:end) = 0;
                case 2 % Itti
                    w(:, 1:13) = 0;
                    w(:, 17:end) = 0;
                case 3 % Color
                    w(1:16) = 0;
                    w(:, 28:end) = 0;
                case 4 % Torralba
                    w(:, 1:27) = 0;
                    w(:, 29:end) = 0;
                case 5 % Horizon
                    w(:, 1:28) = 0;
                    w(:, 30:end) = 0;
                case 6 % Object
                    w(:, 1:29) = 0;
                    w(:, 32:end) = 0;
                case 7 % Center
                    w(:, 1:32) = 0;
                case 8 % All
                    % No modification needed
            end
        
            saliencyMap = (features * w(1:end-1)') + w(end);
            saliencyMap = (saliencyMap - min(saliencyMap)) / (max(saliencyMap) - min(saliencyMap));
            saliencyMap = reshape(saliencyMap, resizedDims);
            saliencyMap = imresize(saliencyMap, [imgWidth, imgHeight]);
        
            [file, type, subjectName];
            X = eyeData(:, 1);
            Y = eyeData(:, 2);
            X(isnan(X)) = [];
            Y(isnan(Y)) = [];
            origImgSize = size(saliencyMap);
            X1 = X(1:floor(end/2));
            X2 = X(floor(end/2)+1:end);
            Y1 = Y(1:floor(end/2));
            Y2 = Y(floor(end/2)+1:end);
            ROCs(subjectName, file, type, 1) = rocScoreSaliencyVsFixations(saliencyMap, X1, Y1, origImgSize);
            ROCs(subjectName, file, type, 2) = rocScoreSaliencyVsFixations(saliencyMap, X2, Y2, origImgSize);
        end

        if fileCounter == 2
            break;
        end
        fileCounter = fileCounter + 1;
    end
end

%% Histogram of ROCs
for featureNum = 1:8
    figure;
    histogram(ROCs(:,:,featureNum,1), 20, 'FaceColor','red');
    hold on;
    histogram(ROCs(:,:,8,1), 20, 'FaceColor','blue');
    title('Comparing ROCs - First Half');
    xlabel('ROC');
    ylabel('Count');
    legendinf = "Feature " + num2str(featureNum);
    legend(legendinf, 'All Features');
    set(gcf, 'PaperPositionMode', 'auto')
    print("D:/Sharif/Term2/Advanced Neuroscience/HW/8/images/saliency/part3_comparingROCs_" + num2str(featureNum) , '-dpng', '-r0')
end

% Mean of ROCs
figure;
meanROCs = squeeze(mean(ROCs, [1 2]));
errorROCs = squeeze(std(ROCs, 0, [1 2])) / sqrt(15 * size(files, 2));

errorbar(1:8, meanROCs(:, 1), errorROCs(:, 1), 'red');
hold on;
errorbar(1:8, meanROCs(:, 2), errorROCs(:, 2), 'blue');
grid on;
grid minor;
legend('Bottom-up', 'Top-down');
title('Mean of ROCs');
ylabel('ROC');
set(gca, 'xtick', 1:8, 'xticklabel', {'Subband', 'Itti', 'Color', 'Torralba', 'Horizon', 'Object', 'DistToCenter', 'All Features'});
xtickangle(90);
set(gcf, 'PaperPositionMode', 'auto')
print("D:/Sharif/Term2/Advanced Neuroscience/HW/8/images/saliency/part3_meanOfROCs" , '-dpng', '-r0')
