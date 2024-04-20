clc;
clear;
close all;

% Frequency range and number of wavelet cycles
frequencies = 4:60;
n_cycles = 4;

% Define the file path template
basePath = 'S001\';
runs = {'R03', 'R04', 'R05', 'R06', 'R07', 'R08', 'R09', 'R10', 'R11', 'R12', 'R13', 'R14'};

% Selected electrodes based on relevance to motor and imagery tasks
electrodes = {'C3__', 'C4__', 'Fcz_', 'Cz__', 'P3__', 'P4__'};
nElectrodes = length(electrodes);

% Assume runs R03, R05, R07, R09, R11, R13 are motor tasks
% and runs R04, R06, R08, R10, R12, R14 are imagery tasks
motorRuns = {'R03', 'R05', 'R07', 'R09', 'R11', 'R13'};
imageryRuns = {'R04', 'R06', 'R08', 'R10', 'R12', 'R14'};

% Initialize storage for PLV results by task type
PLVResultsMotor = struct();
PLVResultsImagery = struct();

% Modify your existing PLV calculation loop:
for r = 1:length(runs)
    filePath = [basePath 'S001' runs{r} '.edf'];
    [data, header] = edfread(filePath);
    for i = 1:nElectrodes
        for j = i+1:nElectrodes
            electrode1 = electrodes{i};
            electrode2 = electrodes{j};
            channelData1 = double(cell2mat(data{:, electrode1}));
            channelData2 = double(cell2mat(data{:, electrode2}));
            plvResults = calculatePLV(channelData1, channelData2, 160, frequencies, n_cycles);
            
            % Categorize by task type
            if ismember(runs{r}, motorRuns)
                if ~isfield(PLVResultsMotor, electrode1) || ~isfield(PLVResultsMotor.(electrode1), electrode2)
                    PLVResultsMotor.(electrode1).(electrode2) = [];
                end
                PLVResultsMotor.(electrode1).(electrode2) = [PLVResultsMotor.(electrode1).(electrode2), plvResults];
            elseif ismember(runs{r}, imageryRuns)
                if ~isfield(PLVResultsImagery, electrode1) || ~isfield(PLVResultsImagery.(electrode1), electrode2)
                    PLVResultsImagery.(electrode1).(electrode2) = [];
                end
                PLVResultsImagery.(electrode1).(electrode2) = [PLVResultsImagery.(electrode1).(electrode2), plvResults];
            end
        end
    end
end

% Preallocate matrix to store p-values
pValues = struct();

for i = 1:nElectrodes
    for j = i+1:nElectrodes
        electrode1 = electrodes{i};
        electrode2 = electrodes{j};
        if isfield(PLVResultsMotor, electrode1) && isfield(PLVResultsMotor.(electrode1), electrode2)
            motorData = PLVResultsMotor.(electrode1).(electrode2);
            imageryData = PLVResultsImagery.(electrode1).(electrode2);
            p = zeros(length(frequencies), 1);
            for k = 1:length(frequencies)
                % Perform t-test at each frequency
                [~, p(k)] = ttest2(motorData(k, :), imageryData(k, :));
            end
            pValues.(electrode1).(electrode2) = p;
        end
    end
end

% Visualization of PLV Results for Multiple Electrode Pairs with Significance
figure;
pairIndex = 1;
totalPairs = nElectrodes * (nElectrodes - 1) / 2;  % Calculate total pairs for subplot arrangement
for i = 1:nElectrodes
    for j = i+1:nElectrodes
        electrode1 = electrodes{i};
        electrode2 = electrodes{j};
        if isfield(pValues, electrode1) && isfield(pValues.(electrode1), electrode2)
            subplot(ceil(totalPairs/3), 3, pairIndex);
            motorPlvValues = PLVResultsMotor.(electrode1).(electrode2);
            imageryPlvValues = PLVResultsImagery.(electrode1).(electrode2);
            plot(frequencies, motorPlvValues, 'b-', 'LineWidth', 2);
            hold on;
            plot(frequencies, imageryPlvValues, 'r-', 'LineWidth', 2);
            % Highlight significant frequencies
            sigFreqs = frequencies(pValues.(electrode1).(electrode2) < 0.05);
            scatter(sigFreqs, motorPlvValues(pValues.(electrode1).(electrode2) < 0.05), 'ko');
            scatter(sigFreqs, imageryPlvValues(pValues.(electrode1).(electrode2) < 0.05), 'ko');
            title(['PLV between ' electrode1 ' and ' electrode2]);
            xlabel('Frequency (Hz)');
            ylabel('Phase Locking Value');
            legend('Motor', 'Imagery', 'Significant');
            grid on;
            pairIndex = pairIndex + 1;
        end
    end
end
sgtitle('Phase Locking Value across Different Electrode Pairs with Significance');
set(gcf, 'Position', [100, 100, 1200, 800]);  % Resize figure to fit all subplots

% Function for PLV calculation
function plvResult = calculatePLV(signal1, signal2, fs, frequencies, n_cycles)
    time = -0.5:1/fs:0.5;
    plvResult = zeros(length(frequencies), 1);
    for f_idx = 1:length(frequencies)
        f = frequencies(f_idx);
        wavelet1 = waveletTransform(signal1, fs, f, n_cycles);
        wavelet2 = waveletTransform(signal2, fs, f, n_cycles);
        phaseDiff = angle(wavelet1) - angle(wavelet2);
        plvResult(f_idx) = abs(mean(exp(1i * phaseDiff)));
    end
end

% Function for wavelet transformation
function waveletResult = waveletTransform(signal, fs, frequency, n_cycles)
    time = -0.5:1/fs:0.5;
    s = n_cycles / (2 * pi * frequency);
    A = 1 / sqrt(s * sqrt(pi));
    wavelet = A * exp(-(time.^2) / (2 * s^2)) .* exp(1i * 2 * pi * frequency * time);
    waveletResult = conv(signal, wavelet, 'same');
end
