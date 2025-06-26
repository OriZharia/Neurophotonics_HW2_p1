
%% Neurophotonics Homework Assignment 2

%% 1. Initialize
% Clear workspace, command window, and close all figures
% Add auxiliary code path

clc;
clear all;
close all;
addpath('./auxiliary_code/');
windowSize = 7;

%% 2. Select Recording Directory
% Prompt user to select a folder containing the recordings

recordingDir = uigetdir('', 'Select a folder');
% Check if a folder was selected
if recordingDir == 0
    error('No folder was selected. Exiting...');
end

%% 3. Discover and Parse Folders
% Get a list of subfolders in the selected directory
% from the folder name, we will extract the parameters and video recordings

subfolders = dir(recordingDir);
subfolders = subfolders([subfolders.isdir]);
% Filter out the current and parent directory entries
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));

% Create a dynamic list of video recordings and their metadata
recordData = {};
videoRecordings = {};
divRec = false;
for i = 1:length(subfolders)
    folderName = subfolders(i).name;
    fullPath = fullfile(recordingDir, folderName);
    fprintf('\nAnalyzing: %s\n', folderName);
    % Check if the folder name contains all required parameters
    if ~contains(folderName, 'Gain') || ~contains(folderName, 'expT') || ~contains(folderName, 'BL') || ~contains(folderName, 'FR')
        warning('Folder "%s" does not contain all required parameters.', folderName);
        continue;
    end
    % Check if the folder contains .tiff files
    tiffFiles = dir(fullfile(fullPath, '*.tiff'));
    if isempty(tiffFiles)
        warning('No .tiff files found in folder: %s', fullPath);
        continue;
    end

    try
        % Read recordings 
        recordings = ReadRecord(fullPath);
        % Check if recordings are empty
        if isempty(recordings)
            warning('No recordings found in folder: %s', fullPath);
            continue;
        end
        % add recordings to the list
        videoRecordings{end + 1} = recordings;
        fprintf('Loaded image successfully from: %s\n', folderName);
        % get the recording info
        info = GetRecordInfo(fullPath);
        % Check if the recordings are divisible by 16
        if divRec || all(mod(recordings(:), 16) == 0)
            recordings = recordings ./ 16;
            divRec = true;
            % info.numBits = 12; % Adjust bit depth if necessary
        end
        % Store the recording metadata
        recordData{end+1} = struct( ...
            'BlackLevel', info.name.BL, ...
            'FrameRate', info.name.FR, ...
            'Gain_dB', info.name.Gain, ...
            'ExposureTime_ms', info.name.expT, ...
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
% Check if dark video index was found
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
% Check if video index was found
if exist('videoIndex', 'var') == 0
    error('No video found in the selected recordings.');
end

%     dark_size = size(videoRecordings{darkVideoIndex}, 3);


%% 4. Define ROI
% Prompt user to define a Region of Interest (ROI) for the video
% If an ROI file already exists, prompt the user to load it or define a new one
% If the ROI is not defined, we will use the first frame of the video to define it


% Check if ROI is already defined
recordingDir= recordData{videoIndex}.FullPath;
roiFile = fullfile(recordingDir, 'ROI_Mask.mat');
if exist(roiFile, 'file')
    % If ROI file exists, prompt user to load it
    fprintf('ROI file found, do you want to load it? (y/n): ');
    userInput = input('', 's');
    % If user chooses not to load, define a new ROI
    if lower(userInput) ~= 'y'
        try
            [mask, roi, figRoi] = ROIMask(videoRecordings{videoIndex}(:,:,1));
            % Save the new ROI to file
            save(roiFile, 'mask', 'roi');
            % save ROi figure as fig
            savefig(figRoi, fullfile(recordingDir, 'ROI_Figure.fig'));
            close(figRoi); % Close the figure after saving

        catch ME
            error('Error defining ROI: %s', ME.message);
        end
    % Load the ROI from the file
    else
        try
            load(roiFile, 'mask', 'roi');
            fprintf('Loaded ROI from file: %s\n', roiFile);
        % If loading fails, prompt user to define a new ROI
        catch ME
            warning('Error loading ROI from file, will define a new ROI instead');
            try
                [mask, roi,figRoi] = ROIMask(videoRecordings{videoIndex}(:,:,1));
                % Save the new ROI to file
                save(roiFile, 'mask', 'roi');
                 % save ROi figure as fig
                savefig(figRoi, fullfile(recordingDir, 'ROI_Figure.fig'));
                close(figRoi); % Close the figure after saving
            catch ME
                error('Error defining ROI: %s', ME.message);
            end
        end
    end
% No ROI file found, prompt user to define a new ROI
else
    fprintf('No ROI file found. Please define a new ROI.\n');
    % Load the first frame of the video to define ROI
    try
        [mask, roi] = ROIMask(videoRecordings{videoIndex}(:,:,1));
        % Save the new ROI to file
        save(roiFile, 'mask', 'roi');
    catch ME
        error('Error defining ROI: %s', ME.message);
    end
end

%% 5. Compute Background Noise (Read Noise)
% Check if read-noise data already exists
% If the read-noise data is already saved, load it
% If not, compute them from the dark video recordings
% Use the dark video index to get the correct recording directory

recordingDir= recordData{darkVideoIndex}.FullPath;
% Check read-noise data already exists
readNoiseFile = fullfile(recordingDir, 'SCOS_Read_Noise.mat');
if exist(readNoiseFile, 'file')
    try
        load(readNoiseFile, 'backgroundImg', 'var_r', 'darkVarPerWindow');
        fprintf('Loaded pixel non-uniformity data from file: %s\n', readNoiseFile);
    catch
        warning('Error loading pixel non-uniformity data from file: %s, will compute them instead.', readNoiseFile);
        backgroundImg = mean(videoRecordings{darkVideoIndex}, 3);
        var_r = std(videoRecordings{darkVideoIndex}(:,:,:), 0, 3).^2;
        darkVarPerWindow = imboxfilt(var_r, windowSize);
        save(readNoiseFile, 'backgroundImg', 'var_r', 'darkVarPerWindow');
        fprintf('Calculated and saved pixel non-uniformity data to file: %s\n', readNoiseFile);
    end
else
    backgroundImg = mean(videoRecordings{darkVideoIndex}, 3);
    var_r = std(videoRecordings{darkVideoIndex}, 0, 3).^2;
    darkVarPerWindow = imboxfilt(var_r, windowSize);
    save(readNoiseFile, 'backgroundImg', 'var_r', 'darkVarPerWindow');
    fprintf('Calculated and saved pixel non-uniformity data to file: %s\n', readNoiseFile);
end


%% 6. Compute Pixel Non-Uniformity
% Check if pixel non-uniformity data already exists
% If the pixel non-uniformity data is already saved, load it
% If not, compute them from the video recordings
% Use the video index to get the correct recording directory

recordingDir= recordData{videoIndex}.FullPath;

% Check pixel non-uniformity data are already saved
PixelNonUniformityFile = fullfile(recordingDir, 'SCOS_Pixel_NonUniformity.mat');
% If the file exists, load the background mean and variance
if exist(PixelNonUniformityFile, 'file')
    try
        % Load the background mean and variance from the file
        load(PixelNonUniformityFile, 'mean_Isp', 'var_sp');
        fprintf('Loaded pixel non uniformity from file: %s\n', PixelNonUniformityFile);
    % If the file does not exist or loading fails, compute them
    catch ME
        warning('Error loading pixel non uniformity from file: %s, will compute them instead.');
        vid_size = size(videoRecordings{videoIndex}, 3);
        if vid_size < 500
            error('Not enough frames in the video. At least 500 frames are required.');
        end
        mean_Isp = mean(videoRecordings{videoIndex}, 3) - recordData{videoIndex}.BlackLevel;
        var_sp = stdfilt(mean_Isp, true(windowSize)).^2;
        save(PixelNonUniformityFile, 'mean_Isp', 'var_sp');
        fprintf('Calculated and saved pixel non uniformity to file: %s\n', PixelNonUniformityFile);
    end 
% there is no background mean and variance file, compute them  
else
    mean_Isp = mean(videoRecordings{videoIndex}, 3) - recordData{videoIndex}.BlackLevel;
    var_sp = stdfilt(mean_Isp, true(windowSize)).^2;
    save(PixelNonUniformityFile, 'mean_Isp', 'var_sp');
    fprintf('Calculated and saved pixel non uniformity to file: %s\n', PixelNonUniformityFile);
end

figure;
imagesc(mean_Isp); colorbar;
title('Mean Background Image');

figure;
imagesc(var_sp); colorbar;
title('Background Variance Image');


%% 7. Estimate System Gain
% Check if gain calculation data already exists
% If the gain calculation data is already saved, load it
% If not, compute them from the video recordings
% Use the video index to get the correct recording directory

recordingDir= recordData{videoIndex}.FullPath;
% Check if gain calculation data already exists
gainFile = fullfile(recordingDir, 'SCOS_Gain_Calculation.mat');
gainFigFile = fullfile(recordingDir, 'SCOS_Gain_Calculation_Fig.png');

% If the gain calculation file exists, load it
if exist(gainFile, 'file')
    try
        load(gainFile, 'gainCalc', 'gainTheoretical');
        fprintf('Loaded gain calculation data from file: %s\n', gainFile);
    % If loading fails, compute them
    catch
        warning('Error loading gain calculation data from file: %s, will compute them instead.', gainFile);
        meanDark = backgroundImg;
        efficiency = 10500; % Example efficiency value, adjust as needed
        [gainCalc, gainTheoretical, gainFig] = GainCalc(videoRecordings{videoIndex}, meanDark, mask, efficiency, recordData{videoIndex}.Gain_dB, recordData{videoIndex}.Bits);
        save(gainFile, 'gainCalc', 'gainTheoretical');
        % Save the gain figure
        saveas(gainFig, gainFigFile);
        close(gainFig); % Close the figure after saving
        fprintf('Calculated and saved gain calculation data to file: %s\n', gainFile);
    end
% If the gain calculation file does not exist, compute them
else
    meanDark = backgroundImg;
    efficiency = 10500; % Example efficiency value, adjust as needed
    [gainCalc, gainTheoretical, gainFig] = GainCalc(videoRecordings{videoIndex}, meanDark, mask, efficiency, recordData{videoIndex}.Gain_dB, recordData{videoIndex}.Bits);
    save(gainFile, 'gainCalc', 'gainTheoretical');
    % Save the gain figure
    saveas(gainFig, gainFigFile);
    close(gainFig); % Close the figure after saving
    fprintf('Calculated and saved gain calculation data to file: %s\n', gainFile);
end



%% 8. Compute K² per Frame

% Initialize variables for K² calculation
numFrames = size(videoRecordings{videoIndex}, 3);
K2_raw = zeros(1, numFrames);
K2_corrected = zeros(1, numFrames);
BFi = zeros(1, numFrames);
% Show progress with a waitbar
w = waitbar(0, 'Calculating K² values...');

for frame = 1:numFrames
    
    im = videoRecordings{videoIndex}(:,:,frame);
    meanIm = mean(im(mask) - backgroundImg(mask))-(recordData{videoIndex}.BlackLevel - recordData{darkVideoIndex}.BlackLevel);
    Var = stdfilt(im, true(windowSize)).^2;

    K2_raw(frame) = mean(Var(mask)) / (meanIm^2);
    K2_corrected(frame) = mean(Var(mask) - gainCalc * meanIm - darkVarPerWindow(mask) - var_sp(mask) - 1/12) / (meanIm^2);
    BFi(frame) = 1 / K2_corrected(frame);
    if mod(frame, 100) == 0
        waitbar(frame / numFrames, w, sprintf('Calculating K² values... Frame %d of %d', frame, numFrames));
    end
end
waitbar(1, w, 'Calculation complete!'); % Update waitbar to show completion
close(w); % Close the waitbar

%% 9. Plot K² Over Time
% Create a time vector based on the frame rate of the video
timeVector = (1:numFrames) / recordData{videoIndex}.FrameRate;
% Get the recording directory for saving plots
recordingDir = recordData{videoIndex}.FullPath;
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
savefig(gcf,fullfile(recordingDir, 'SCOS_K2_Over_Time.fig'))

% Plot BFi in a separate plot
figure;
plot(timeVector, BFi, 'g-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('BFi');
title('BFi Over Time');
grid on;
saveas(gcf, fullfile(recordingDir, 'SCOS_BFi_Over_Time.png'));
savefig(gcf, fullfile(recordingDir, 'SCOS_BFi_Over_Time.fig'));


%% 10. Save Results and figures
% Save all results in a structured format
resultsFile = fullfile(recordingDir, 'SCOS_Results.mat');
save(resultsFile,'timeVector', 'K2_raw', 'K2_corrected', 'BFi');
fprintf('Results saved to: %s\n', resultsFile);
% close all figures
close all;

