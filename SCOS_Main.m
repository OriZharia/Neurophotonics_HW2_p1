%% Recording Integrity Check 
%%
% Opens a dialog to select the root recordings folder.
% User must manually select the "recordings" folder (one level up).
% Scans its subfolders.
% Looks for tiff files in each subfolder.
% Checks if the first image can be read.
% Prints warnings if there are problems.
% Checks if all pixels are divisible by 16 (MSB check), as required in the assignment.

%%
%%%% Should we also check the contents of each folder and do excessive validation and add defaults, etc.?
%%%% Not sure if it's necessary

%%
%% Recording Integrity Check 
% clear;
clc;
clear all;
% clear figure;
close all;
addpath('./auxiliary_code/');
%% check validity of recordings


%% Select Recording Directory
recordingDir = uigetdir('C:\Users\ayele\Desktop\תואר שני\נוירופוטוניקה\מטלה 2', ...
                        'Select a folder');
if recordingDir == 0
    error('No folder was selected. Exiting...');
end

%% Discover Subfolders
subfolders = dir(recordingDir);
subfolders = subfolders([subfolders.isdir]);
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));




%% Loop Through Folders: Check .tiff and Extract Metadata
recordData = cell(length(subfolders), 1);
videoRecordings = cell(length(subfolders), 1);
divRec=false;
for i = 1:length(subfolders)
    folderName = subfolders(i).name;
    fullPath = fullfile(recordingDir, folderName);
    fprintf('\nAnalyzing: %s\n', folderName);

    % Check for required keywords in folder name
    if ~contains(folderName, 'Gain') || ...
       ~contains(folderName, 'expT') || ...
       ~contains(folderName, 'BL') || ...
         ~contains(folderName, 'FR')
        warning('Folder "%s" does not contain all required parameters.', folderName);
        continue;
    end

    % Check if .tiff files exist
    tiffFiles = dir(fullfile(fullPath, '*.tiff'));
    if isempty(tiffFiles)
        warning('No .tiff files found in folder: %s', fullPath);
        continue;
    end

    % Try to load image and metadata
    try
        recordings = ReadRecord(fullPath);
        videoRecordings{i} = recordings;
        fprintf('Loaded image successfully from: %s\n', folderName);
        info = GetRecordInfo(fullPath);
        if divRec || all(mod(recordings(1:400),16) == 0)
            rec = rec/16; % because Basler camera for some reason uses last 12 bits instead of first
            divRec=true;
        end
        

        
       recordData{i} = struct( ...
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

%% Get dark video index 
% Remove empty entries from videoRecordings and recordData
for vids = 1:length(videoRecordings)
    if isempty(videoRecordings{vids})
        videoRecordings(vids) = [];
        recordData(vids) = [];
    end
end

% Find the first dark video and the first video (non-dark)
darkVideoIndex = [];
videoIndex = [];
for i = 1:length(recordData)
    if contains(recordData{i}.FolderName,'dark')
        darkVideoIndex(end+1) = i;  % Store index of dark video
    else 
        videoIndex(end+1) = i;
    end
end
% Check if dark video and video exist
if isempty(darkVideoIndex)
    error('No dark video found in the recordings.');
end
if isempty(videoIndex)
    error('No video found in the recordings.');
end
% check if there are more than one dark video or video
if length(darkVideoIndex) > 1
    warning('Multiple dark videos found. Using the first one.');
    darkVideoIndex = darkVideoIndex(1);
end
if length(videoIndex) > 1
    warning('Multiple videos found. Using the first one.');
    videoIndex = videoIndex(1);
end


%% get mask
[mask, roi]= ROIMask(videoRecordings{videoIndex}(:,:,1));


%% שלב 5 – חישוב רקע: σ_bg ו- <I_bg>



% חישוב סטטיסטיקות
mean_Ibg = mean(videoRecordings{darkVideoIndex},3);
sigma_bg = std(mean_Ibg(mask));

% הדפסה
% fprintf('<I_bg> = %.3f\n', mean_Ibg);
fprintf('σ_bg = %.3f\n', sigma_bg);


%% Calculate the Gain
meanDark=mean_Ibg;
efficiency=10500;
recordData{videoIndex}.Bits=12;
[gainCalc, gainTheotycal ]= GainCalc(videoRecordings{videoIndex}, meanDark,mask, efficiency, recordData{videoIndex}.Gain_dB, recordData{videoIndex}.Bits);