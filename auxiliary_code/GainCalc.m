
function [gainCalc, gainTheotycal ]= GainCalc(recVid, meanDark,mask, efficiency, gainVid, numBits)
% GainCalc calculates the gain of a system based on amplitude, frequency, and phase.
% It uses a linear fit to the variance and mean of the video frames.
% Inputs:
%   recVid: Video record (cell array of video frames)
%   meanDark: Mean dark frame (2D matrix)
%   mask: ROI mask to apply on the frames
%   efficiency: Efficiency of the camera system
%   gainVid: Gain value from the video metadata (in dB)
%   numBits: Number of bits of the camera system

% Outputs:
%   gainCalc: Calculated gain of the camera system
%   gainTheotycal: Theoretical gain based on camera system parameters



%% Get the video length
videoLength = size(recVid, 3);

%% Calcualre Variance and mean
meanFrame = zeros(videoLength, 1);
varFrame = zeros(videoLength, 1);
for i = 1:(videoLength)
    rec= recVid(:,:,i);
    meanFrame(i) = mean(rec(mask)-meanDark(mask), "all");
    varFrame(i) = var(rec(mask)- meanDark(mask));
end


%% Create a linear fit
% Var(I)=G/e * Mean(I) + noise^2
p = polyfit(meanFrame, varFrame, 1);
%% Calculate gain by using the slope of the linear fit
% Gain is the slope of the linear fit divided by the efficiency
gainBase= p(1);
% gainBase is the gain set by the camera system 

gainIn=  10^(gainVid/20);
% Calculate the gain of the camera system
gainCalc= gainBase;

gainTheotycal = 2^numBits / efficiency * gainIn;

%% Create a figure to visualize the fit with the calculated gain
figure;
scatter(meanFrame, varFrame, 'filled');
hold on;
x = linspace(min(meanFrame), max(meanFrame), 10);
y = polyval(p, x);
plot(x, y, 'r-', 'LineWidth', 2);
xlabel('Mean Intensity');
ylabel('Variance');
title('Variance vs Mean Intensity with Linear Fit');
grid on;
% Display the calculated gain inside the plot
text(0.05 * (max(meanFrame)-min(meanFrame)) + min(meanFrame), ...
    0.95 * (max(varFrame)-min(varFrame)) + min(varFrame), ...
    {['Calculated Gain: ', num2str(gainCalc, '%.3g')], ...
     ['Theoretical Gain: ', num2str(gainTheotycal, '%.3g')]}, ...
    'FontSize', 10, 'VerticalAlignment', 'top', 'BackgroundColor', 'w');

%% save the figure
saveas(gcf, 'GainCalculation.png');
% Close the figure
close(gcf);
end



% function G= GainCalc(filepaths, meanDark,BLDark, efficiency, numBits)
% % GainCalc calculates the gain of a system based on amplitude, frequency, and phase.
% % It uses a linear fit to the variance and mean of the video frames.
% % Inputs:
% %   filepaths: Cell array of file paths to video files
% %   meanDark: Mean dark frame
% %   BLDark: Black level of the dark
% %   efficiency: Efficiency of the camera system
% %   numBits: Number of bits of the camera system
% % Outputs:
% %   G: Calculated gain of the camera system
% %   gainTheotycal: Theoretical gain based on camera system parameters
% % Initialize variables
% meanVec = zeros(size(filepaths, 2), 1);
% varVec = zeros(size(filepaths, 2), 1);
% varErrorVec = zeros(size(filepaths, 2), 1);
% 
% % Loop through each video file
% for i = 1:length(filepaths)
%     % Read the video file
%     recVid = ReadRecord(filepaths{i});
% 
%     % Get the number of frames in the video
%     numFrames = size(recVid, 3);
%     info=GetRecordInfo(filepaths{i});
%     BLVid= info.name.BlackLevel;
%     gainVid=info.name.Gain;
%     gainVidLin=10^(gainVid/20);
%     meanFrame = zeros(numFrames, 1);
%     varFrame = zeros(numFrames, 1);
% 
% 
%     % Loop through each frame in the video
%     for j = 1:numFrames
%         frame = recVid(:, :, j);
%         % Calculate mean and variance for the current frame
%         meanFrame(j) = (mean(frame(:)-meanDark(:)) - (BLVid - BLDark))/gainVidLin;
%         varFrame(j) = var(frame(:)-meanDark(:));
%     end
%     % Store the mean and variance for the current video
%     meanVec(i) = mean(meanFrame);
%     varVec(i) = mean(varFrame);
%     varErrorVec(i) = std(varFrame) / sqrt(numFrames); % Standard error of the variance
% 
% end
% %% Account for the gain
% % Create a linear fit
% p = polyfit(meanVec, varVec, 1);
% % Calculate gain by using the slope of the linear fit
% gainBase = p(1);
% % Gain is the slope of the linear fit divided by the efficiency
% G = gainBase;
% % Theoretical gain based on camera system parameters
% gainTheotycal = 2^numBits / efficiency;
% % Create a figure to visualize the fit with the calculated gain
% figure;
% scatter(meanVec, varVec, 'filled');
% hold on;
% errorbar(meanVec, varVec, varErrorVec, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'CapSize', 8);
% x = linspace(min(meanVec), max(varVec), 100);
% y = polyval(p, x);
% plot(x, y, 'r-', 'LineWidth', 2);
% xlabel('Mean Intensity');
% ylabel('Variance');
% title('Variance vs Mean Intensity with Linear Fit');
% grid on;
% hold off;
% % Display the calculated gain
% disp(['Calculated Gain: ', num2str(G)]);
% disp(['Theoretical Gain: ', num2str(gainTheotycal)]);
% 
% end
% 