%% Preprocess
% Preprocesses TMS-EEG data from https://download.fieldtriptoolbox.org/tutorial/tms/sp/
%   to be more convenient to load into an R environment and
%   more useful outside of the specific structure perscribed by the Field-Trip MATLAB toolbox

%% Read data
addpath("src")
data_dir = fullfile("data");
addpath(data_dir);
dataset = load("data_tms_clean.mat");

disp("Before:")
dataset = dataset.data_tms_clean

%% Create time-locked averages
% https://www.fieldtriptoolbox.org/tutorial/tms/tms-eeg/#time-locked-averaging

addpath(genpath("lib/fieldtrip-20261111"))

cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.05 -0.001];
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq = 35;

% Find all trials corresponding to the relax condition
cfg.trials = find(dataset.trialinfo==1);
relax_avg = ft_timelockanalysis(cfg, dataset);

% Find all trials corresponding to the contract condition
cfg.trials = find(dataset.trialinfo==3);
contract_avg = ft_timelockanalysis(cfg, dataset);

% Calculate the difference ERP wave
cfg = [];
cfg.operation = 'subtract'; % Operation to apply
cfg.parameter = 'avg'; % The field in the data structure to which to apply the operation
difference_avg = ft_math(cfg, contract_avg, relax_avg);

% Plot TEPs of both conditions
cfg = [];
cfg.layout = 'easycapM10'; % Specifying this allows you to produce topographical plots of your data
cfg.channel = '17';
cfg.xlim = [-0.1 0.6];
ft_singleplotER(cfg, relax_avg, contract_avg, difference_avg);
ylabel('Amplitude (uV)');
xlabel('time (s)');
title('Relax vs Contract');
legend({'relax' 'contract' 'contract-relax'});

%% Strip dataset to only necessary information

% Remove unnecessary fields
redundant_fields = ["hdr", "cfg", "fsample"];

for f = redundant_fields
    if ~isfield(dataset, f)
        continue
    end
    dataset = rmfield(dataset, f);
end

% De-duplicate time vectors
dataset.time = dataset.time{1};

% Convert channel labels to integers
dataset.label = int32(str2double(string(dataset.label)));

% Concatenate trials into multidimensional array
dataset.data = cat(3, dataset.trial{:});
dataset = rmfield(dataset, "trial");

% Add information
dataset.trialseq = int32(1:size(dataset.data, 3));

disp("After:")
disp(dataset)

%% Save dataset
% Save using HDF5 format (v7.3)
save(fullfile(data_dir, "data_tms_clean_for_R.mat"), "-struct", "dataset", "-v7.3")