function [mask, roi] = ROIMask(fig)
    % Function to create a circular ROI mask on a 2D image
    % Input:
    %   fig - 2D image matrix
    % Outputs:
    %   mask - binary mask with 1s inside the circle, 0s outside
    %   roi - circle ROI object (contains center and radius info)

    figure;
    imagesc(fig);
    axis image;
    colormap gray;
    title('Draw a Circular ROI');
    
    roi = drawcircle('Color', 'r');  % Interactive circle ROI
    
    % Create coordinate grid
    [X, Y] = meshgrid(1:size(fig,2), 1:size(fig,1));
    
    % Compute distance of each pixel to the ROI center
    dist_from_center = sqrt((X - roi.Center(1)).^2 + (Y - roi.Center(2)).^2);
    
    % Create mask: 1 inside circle, 0 outside
    mask = dist_from_center <= roi.Radius;
end