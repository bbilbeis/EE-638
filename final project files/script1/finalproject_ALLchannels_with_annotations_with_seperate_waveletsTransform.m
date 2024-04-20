clc;
clear;
close all;

filePath = 'S001\S001R14.edf'; % path of eeg data for subject 1 run 14

% Read the EDF file
[data, header] = edfread(filePath);

% Initialization for annotations and annotation codes
annotations = header;
onsets = seconds(annotations.Onset);
types = annotations.Annotations;
durations = seconds(annotations.Duration);

fs = 160; % Sampling Frequency
frequencies = 4:60; % 4 to 60 Hz (wavelet analysis)
n_cycles = 4; % wavelet cycles

% Initialize results container
numChannels = length(data.Properties.VariableNames);
convolutionResults = cell(1, numChannels);

% Process each channel
for i = 1:numChannels
    channelData = data{:, i};

    if iscell(channelData) % convert channelData to numeric if not already
        channelData = cell2mat(channelData);
    elseif isnumeric(channelData)
    else
        error('Channel data is not numeric.');
    end

    % wavelet convolution for data
    convolutionResults{i} = waveletConvolution(channelData, fs, frequencies, n_cycles); 
    
    % Plot Time-Frequency representation
    figure;
    t = linspace(0, length(channelData) / fs, length(channelData));
    imagesc(t, frequencies, abs(convolutionResults{i}));
    axis xy;
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
    title(['Time-Frequency Representation of Channel ', data.Properties.VariableNames{i}]);
    colormap('jet');
    colorbar;
end

% plot the combined channel data (post-processing)
figure('Position', [100, 100, 1200, 600]);
hold on;
grid on;
xlabel('Time (sec)');
ylabel('Channels');
title('EEG Data from All Channels');

% vertical offset for each channel
vertical_offsets = linspace(1, numChannels * 15, numChannels);

% Plot each channel
for i = 1:numChannels
    channelData = data{:, i};
    
    if iscell(channelData) % convert channelData to numeric if not already
        channelData = cell2mat(channelData);
    elseif isnumeric(channelData)
        % Data is already numeric
    else
        error('Channel data is not numeric.');
    end

    normalizedChannelData = (channelData - mean(channelData)) / std(channelData);
    t = linspace(0, length(channelData) / 160, length(channelData));
    plot(t, normalizedChannelData + vertical_offsets(i), 'LineWidth', 0.5);
end

% Plot annotations and annotation codes
for j = 1:height(annotations)
    xline(onsets(j), 'Color', 'r', 'LineWidth', 1);
    text(onsets(j), vertical_offsets(end) + 10, types{j}, 'Rotation', 90, 'Color', 'r', 'FontSize', 8);
end

ylim([0, vertical_offsets(end) + 20]); % expand y-axis limits for readability
set(gca, 'ytick', vertical_offsets, 'yticklabel', data.Properties.VariableNames); % show channel names 
hold off;

% Wavelet convolution function
function waveletResult = waveletConvolution(signal, fs, frequencies, n_cycles)
    time = -0.5:1/fs:0.5;
    waveletResult = zeros(length(frequencies), length(signal));
    for f = frequencies
        s = n_cycles / (2 * pi * f);
        A = 1 / sqrt(s * sqrt(pi));
        wavelet = A * exp(-(time.^2) / (2 * s^2)) .* exp(1i * 2 * pi * f * time);
        signalConvolved = conv(signal, wavelet, 'same');
        waveletResult(find(frequencies==f), :) = signalConvolved;
    end
end