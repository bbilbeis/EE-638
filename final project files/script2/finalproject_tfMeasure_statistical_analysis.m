clc;
clear;
close all;

% Initializations
frequencies = 4:60; % Frequency range
n_cycles = 4; % wavelet cycles
basePath = 'S001\';
runs = {'R03', 'R04', 'R05', 'R06', 'R07', 'R08', 'R09', 'R10', 'R11', 'R12', 'R13', 'R14'};
avgMotor = zeros(length(frequencies), 1);
avgImagery = zeros(length(frequencies), 1);
fs = 160; % Sampling frequency

% Loop through each run
for r = 1:length(runs)
    filePath = [basePath 'S001' runs{r} '.edf'];
    [data, ~] = edfread(filePath); % Read the EDF file

    % Process each channel
    for i = 1:length(data.Properties.VariableNames)
        channelData = data{:, i};
        if iscell(channelData)
            channelData = cell2mat(channelData);
        end
        convolutionResult = waveletConvolution(channelData, fs, frequencies, n_cycles);
        if mod(r, 2) == 1 % Real movement tasks
            avgMotor = avgMotor + abs(convolutionResult);
        else % Imagery tasks
            avgImagery = avgImagery + abs(convolutionResult);
        end
    end
end

% Average the results over all runs
avgMotor = avgMotor / (length(runs)/2);
avgImagery = avgImagery / (length(runs)/2);

% Perform t-test to compare motor and imagery tasks
[h, p] = ttest2(avgMotor', avgImagery');

% Visualization
figure;
hold on;
motorLine = plot(frequencies, mean(avgMotor, 2), 'b', 'LineWidth', 2);
imageryLine = plot(frequencies, mean(avgImagery, 2), 'r', 'LineWidth', 2);
significantPoints = scatter(frequencies(p < 0.05), mean(avgMotor(p < 0.05, 2), 2), 'k*');
legend([motorLine, imageryLine, significantPoints], {'Motor', 'Imagery', 'Significant Difference'});
xlabel('Frequency (Hz)');
ylabel('Power');
title('Comparison of Motor and Imagery Tasks');
hold off;

function waveletResult = waveletConvolution(signal, fs, frequencies, n_cycles)
    time = -0.5:1/fs:0.5;
    waveletResult = zeros(length(frequencies), length(signal));
    for f = frequencies
        s = n_cycles / (2 * pi * f);
        A = 1 / sqrt(s * sqrt(pi));
        wavelet = A * exp(-(time.^2) / (2 * s^2)) .* exp(1i * 2 * pi * f * time);
        waveletResult(frequencies == f, :) = conv(signal, wavelet, 'same');
    end
end
