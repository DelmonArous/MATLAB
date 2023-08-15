function [img_dist, xy_closestpeak] = EstimatePeakDistanceMap(bw, ...
    xy_centorids, xy_peak)

% Extract and store centroid coordinates
x_centroids = xy_centorids(:,1);
y_centroids = xy_centorids(:,2);

% Initialization
img_dist = zeros(size(bw,1), size(bw,2));
xy_closestpeak = [];

% Loop over each centroid coordinate
for i = 1:length(x_centroids)

    % Compute the distance from each centroid index i to each peak coordinate
    distances = sqrt((x_centroids(i) - xy_peak(:,1)).^ 2 + ...
        (y_centroids(i) - xy_peak(:, 2)).^2);

    % Find minimum distance (closest peak) to each centroid
    [peakDistance, indexOfMin] = min(distances);

    img_dist(y_centroids(i),x_centroids(i)) = peakDistance;
    xy_closestpeak = [xy_closestpeak; ...
        xy_peak(indexOfMin,1) xy_peak(indexOfMin,2)];

end

end