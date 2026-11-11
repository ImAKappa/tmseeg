%% Preprocess
% Preprocesses TMS-EEG data from https://download.fieldtriptoolbox.org/tutorial/tms/sp/
%   to be more convenient to load into an R environment

%% Read data
addpath("src")
data_dir = fullfile("data");
addpath(data_dir);
dataset = load("data_tms_clean.mat");

disp("Before:")
dataset = dataset.data_tms_clean

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