% Add auxiliary code path
addpath('auxiliary_code\');

% Define all file paths and their labels
filePaths = {
    'C:\Users\Shira\אוניברסיטה\סמסטר ו\מבוא לניורופוטוניקה\HW2\records\4 WhitePaper_Gain24dB_expT0.5ms_BlackLevel0DU', ...
    'C:\Users\Shira\אוניברסיטה\סמסטר ו\מבוא לניורופוטוניקה\HW2\records\1 WithCover_Gain0dB_expT0.021ms_BlackLevel30DU', ...
    'C:\Users\Shira\אוניברסיטה\סמסטר ו\מבוא לניורופוטוניקה\HW2\records\2 WithCover_Gain24dB_expT0.021ms_BlackLevel30DU', ...
    'C:\Users\Shira\אוניברסיטה\סמסטר ו\מבוא לניורופוטוניקה\HW2\records\3 WithCover_Gain24dB_expT10ms_BlackLevel30DU'
};

fileLabels = {
    'With Paper', ...
    'Closed Low Gain, Low Time', ...
    'Closed High Gain, Low Time', ...
    'Closed High Gain, High Time'
};
%constants
capacity=10500;


% Initialize arrays to store results
temporalNoiseVals = zeros(1, length(filePaths));
globalSpatialNoiseVals = zeros(1, length(filePaths));
localSpatialNoisAfterAvgeVals = zeros(1, length(filePaths));
localSpatialNoisBeforeAvgeVals = zeros(1, length(filePaths));
totalNoise=zeros(1, length(filePaths));
temporalNoiseValsSquareDivMean=zeros(1, length(filePaths));
gVal=zeros(1, length(filePaths));


% Analysis loop
for i = 1:length(filePaths)
    filePath = filePaths{i};
    label = fileLabels{i};
    
    % Load record
    rec = ReadRecord(filePath);
    info= GetRecordInfo(filePath);
    Gain=info.name.Gain;
    numBits=info.nBits;
    gVal(i)=(2^numBits)/(capacity)*10^(Gain/20);


    %% Temporal Noise (mean over time)
    TemporalNoise = std(rec, 0, 3);
    temporalNoiseVals(i) = mean(TemporalNoise(:));

    %% Global Spatial Noise (std of time-averaged image)
    meanImage = mean(rec, 3);
    globalSpatialNoiseVals(i) = std(meanImage(:));

    %% Local Spatial Noise (std over local neighborhood after time avg)
    windowSize = 7;
    localStdMap = stdfilt(meanImage, true(windowSize));
    localSpatialNoisAfterAvgeVals(i) = mean(localStdMap(:));
    %% Local Spatial Noise BEFORE time averaging
    windowSize = 7;
    [numRows, numCols, numFrames] = size(rec);
    localStdVals = zeros(1, numFrames);

    for f = 1:numFrames
        frame = rec(:, :, f);
        localStdMap = stdfilt(frame, true(windowSize));
        localStdVals(f) = mean(localStdMap(:));  % Mean of local std in this frame
    end

    localSpatialNoisBeforeAvgeVals(i) = mean(localStdVals);  % Average over all frames
    %% Total Noise
    totalNoise(i)=sqrt(temporalNoiseVals(i)^2+localSpatialNoisAfterAvgeVals(i)^2);
    %% TemporalNoise^2/Mean
    temporalNoiseValsSquareDivMean(i) = mean(TemporalNoise(:).^2 ./ meanImage(:));

end

% Create summary table
NoiseTable = table(fileLabels', temporalNoiseVals', globalSpatialNoiseVals', localSpatialNoisAfterAvgeVals',localSpatialNoisBeforeAvgeVals', ...
    totalNoise',temporalNoiseValsSquareDivMean',gVal','VariableNames', {'File', 'Temporal Noise', 'Global Spatial Noise', 'Local Spatial Noise After Avg', ...
    'Local Spatial Noise Before Avg','Total Noise','Temporal Noise Vals Square Div Mean','G'});
 
% Display table
disp(NoiseTable)

writetable(NoiseTable, 'NoiseResults.csv');

%% plot  imagesc(mean(rec,3)) ) for each file in one figure
% adjust figure size
figure('Position', [100, 100, 1200, 800]);
% Create subplots for each record
for i = 1:length(filePaths)
    rec = ReadRecord(filePaths{i});
    subplot(2, 2, i);
    imagesc(mean(rec, 3)); % Display mean image for each record
    title(fileLabels{i});
    colorbar;
end
% save figure
saveas(gcf, 'MeanImages.png');
% close figure
close all;



