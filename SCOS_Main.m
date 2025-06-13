
%% Neurophotonics Homework Assignment 2

%% 1. Initialize
clc;
clear all;
close all;
addpath('./auxiliary_code/');
windowSize = 7;

%% 2. Select Recording Directory
recordingDir = uigetdir('', 'Select a folder');
if recordingDir == 0
    error('No folder was selected. Exiting...');
end

%% 3. Discover and Parse Folders
subfolders = dir(recordingDir);
subfolders = subfolders([subfolders.isdir]);
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));

% Create a dynamic list of video recordings and their metadata

recordData = {};
videoRecordings = {};
divRec = false;
for i = 1:length(subfolders)
    folderName = subfolders(i).name;
    fullPath = fullfile(recordingDir, folderName);
    fprintf('\nAnalyzing: %s\n', folderName);

    if ~contains(folderName, 'Gain') || ~contains(folderName, 'expT') || ~contains(folderName, 'BL') || ~contains(folderName, 'FR')
        warning('Folder "%s" does not contain all required parameters.', folderName);

        continue;
    end

    tiffFiles = dir(fullfile(fullPath, '*.tiff'));
    if isempty(tiffFiles)
        warning('No .tiff files found in folder: %s', fullPath);
        continue;
    end

    try
        recordings = ReadRecord(fullPath);
        if isempty(recordings)
            warning('No recordings found in folder: %s', fullPath);
            continue;
        end
        videoRecordings{end + 1} = recordings;
        fprintf('Loaded image successfully from: %s\n', folderName);
        info = GetRecordInfo(fullPath);
        if divRec || all(mod(recordings(:), 16) == 0)
            rec = rec / 16;
            divRec = true;
            info.numBits = 12; % Adjust bit depth if necessary
        end

        recordData{end+1} = struct( ...
            'BlackLevel', info.name.BL, ...
            'FrameRate', info.name.FR, ...
            'Gain_dB', info.name.Gain, ...
            'ExposureTime', info.name.expT, ...
            'FileType', info.fileType, ...
            'CameraSN', info.cameraSN, ...
            'Bits', info.nBits, ...
            'FolderName', folderName, ...
            'FullPath', fullPath ...
        );

    catch ME
        error('Error processing folder "%s": %s', folderName, ME.message);
 
    end
end
%% find dark video index
% since this is a cell array, we need to find the index of the first dark video
for i = 1:length(videoRecordings)
    if contains(recordData{i}.FolderName, 'dark')
        darkVideoIndex = i;
        break;
    end
end
if exist('darkVideoIndex', 'var') == 0
    error('No dark video found in the selected recordings.');
end
% find video index
for i = 1:length(videoRecordings)
    if ~contains(recordData{i}.FolderName, 'dark')
        videoIndex = i;
        break;
    end
end
if exist('videoIndex', 'var') == 0
    error('No video found in the selected recordings.');
end



%% 4. Define ROI
try
    [mask, roi] = ROIMask(videoRecordings{videoIndex}(:,:,1));
catch ME
    error('Error defining ROI: %s', ME.message);
end

%% 5. Compute Background Mean and Variance
mean_Ibg = mean(videoRecordings{videoIndex}, 3) - recordData{videoIndex}.BlackLevel;
sigma_bg = stdfilt(mean_Ibg, true(windowSize)).^2;

figure;
imagesc(mean_Ibg); colorbar;
title('Mean Background Image');

figure;
imagesc(sigma_bg); colorbar;
title('Background Variance Image');

%% 6. Compute Pixel Non-Uniformity
dark_size = size(videoRecordings{darkVideoIndex}, 3);
if dark_size < 500
    error('Not enough frames in the dark video. At least 500 frames are required.');
end
backgroundImg = mean(videoRecordings{darkVideoIndex}(:,:,1:500), 3);
std_sp = std(videoRecordings{darkVideoIndex}(:,:,1:500), 0, 3).^2;
darkVarPerWindow = imboxfilt(std_sp, windowSize);

%% 7. Estimate System Gain
meanDark = mean_Ibg;
efficiency = 10500;
[gainCalc, gainTheoretical] = GainCalc(videoRecordings{videoIndex}, meanDark, mask, efficiency, recordData{videoIndex}.Gain_dB, recordData{videoIndex}.Bits);

%% 8. Compute K² per Frame
numFrames = size(videoRecordings{videoIndex}, 3);
K2_raw = zeros(1, numFrames);
K2_corrected = zeros(1, numFrames);
BFi = zeros(1, numFrames);

for frame = 1:numFrames
    im = videoRecordings{videoIndex}(:,:,frame);
    meanIm = mean(im(mask) - backgroundImg(mask)) - (recordData{videoIndex}.BlackLevel - recordData{darkVideoIndex}.BlackLevel);
    Var = stdfilt(im, true(windowSize)).^2;

    K2_raw(frame) = mean(Var(mask)) / (meanIm^2);
    K2_corrected(frame) = mean(Var(mask) - gainCalc * meanIm - darkVarPerWindow(mask) - sigma_bg(mask) - 1/12) / (meanIm^2);
    BFi(frame) = 1 / K2_corrected(frame);
end

%% 9. Plot K² Over Time
timeVector = (1:numFrames) / recordData{videoIndex}.FrameRate;

% Plot K^2 (raw and corrected) in one plot
figure;
subplot(2, 1, 1);
plot(timeVector, K2_raw, 'r-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('K^2');
title('Raw Contrast (K^2) Over Time');
grid on;
subplot(2, 1, 2);
plot(timeVector, K2_corrected, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('K^2');
title('Corrected Contrast (K^2) Over Time');
grid on;
saveas(gcf, fullfile(recordingDir, 'SCOS_K2_Over_Time.png'));

% Plot BFi in a separate plot
figure;
plot(timeVector, BFi, 'g-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('BFi');
title('BFi Over Time');
grid on;
saveas(gcf, fullfile(recordingDir, 'SCOS_BFi_Over_Time.png'));

%% 10. Save Noise Components and Results
noiseFile = fullfile(recordingDir, 'SCOS_Noise_Components.mat');
save(noiseFile, 'mean_Ibg', 'sigma_bg', 'std_sp', 'darkVarPerWindow', 'gainCalc', 'K2_raw', 'K2_corrected', 'BFi');
disp('Noise components and results saved.');
