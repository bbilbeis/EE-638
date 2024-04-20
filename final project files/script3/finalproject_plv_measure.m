clc;
clear;
close all;

frequencies = 4:60; % frequency range
n_cycles = 4; % number of wavelet cycles

% paths and runs
basePath = 'S001\';
runs = {'R03', 'R04', 'R05', 'R06', 'R07', 'R08', 'R09', 'R10', 'R11', 'R12', 'R13', 'R14'};

electrodes = {'C3__', 'C4__', 'Cpz_', 'Cz__', 'Pz__', 'Iz__'}; % electrodes of interest
nElectrodes = length(electrodes);

PLVResults = struct(); % Initialization for PLV results

for r = 1:length(runs)
    filePath = [basePath 'S001' runs{r} '.edf'];
    [data, header] = edfread(filePath);

    for i = 1:nElectrodes % Calculate PLV for each pair of selected electrodes
        for j = i+1:nElectrodes
            electrode1 = electrodes{i};
            electrode2 = electrodes{j};

            % convert to numeric array
            channelData1 = double(cell2mat(data{:, electrode1}));
            channelData2 = double(cell2mat(data{:, electrode2}));

            plvResults = calculatePLV(channelData1, channelData2, 160, frequencies, n_cycles); % plv computation
            PLVResults.(electrode1).(electrode2)(:, r) = plvResults; % store in create PLV results earlier
        end
    end
end

% Average the PLV results over all runs
for i = 1:nElectrodes
    for j = i+1:nElectrodes
        electrode1 = electrodes{i};
        electrode2 = electrodes{j};
        PLVResults.(electrode1).(electrode2) = mean(PLVResults.(electrode1).(electrode2), 2);
    end
end

% plots of PLV Results
figure;
pairIndex = 1;
totalPairs = nElectrodes * (nElectrodes - 1) / 2;  % Calculate total pairs to arrange the subplots properly
for i = 1:nElectrodes
    for j = i+1:nElectrodes
        electrode1 = electrodes{i};
        electrode2 = electrodes{j};
        if isfield(PLVResults, electrode1) && isfield(PLVResults.(electrode1), electrode2)
            subplot(ceil(totalPairs/3), 3, pairIndex);  % Arrange in a grid that approximates square by having three columns
            plot(frequencies, PLVResults.(electrode1).(electrode2), 'LineWidth', 2);
            title(['PLV between ' electrode1 ' and ' electrode2]);
            xlabel('Frequency (Hz)');
            ylabel('Phase Locking Value');
            grid on;
            pairIndex = pairIndex + 1;
        end
    end
end
sgtitle('Phase Locking Value across Different Electrode Pairs');
set(gcf, 'Position', [100, 100, 1200, 800]);  % Resize figure to fit all of the subplots

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

% Function for the wavelet transformation
function waveletResult = waveletTransform(signal, fs, frequency, n_cycles)
    time = -0.5:1/fs:0.5;
    s = n_cycles / (2 * pi * frequency);
    A = 1 / sqrt(s * sqrt(pi));
    wavelet = A * exp(-(time.^2) / (2 * s^2)) .* exp(1i * 2 * pi * frequency * time);
    waveletResult = conv(signal, wavelet, 'same');
end