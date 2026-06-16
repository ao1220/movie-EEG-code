%% LOOCV nested feature selection for SVM

clear; clc;

%% 1. Settings

outputpath = '';
if ~exist(outputpath, 'dir')
    mkdir(outputpath);
end

% Set these folders before running.
inputpath_MCI_m1 = '';
inputpath_MCI_m2 = '';
inputpath_HC_m1  = '';
inputpath_HC_m2  = '';

alpha = 0.05;

% Options: 'uncorrected', or 'candidate_edges_fdr'.
fdrMode = 'candidate_edges_fdr';

useBands = 1:5;

band_name = {'delta','theta','alpha','beta','gamma','fullband'};

feature_names = { ...
    'negative_delta','negative_theta','negative_alpha','negative_beta','negative_gamma', ...
    'positive_delta','positive_theta','positive_alpha','positive_beta','positive_gamma'};

%% 2. Load data

disp('Loading HC movie1...');
[network_HC_m1, files_HC_m1] = load_network_dir(inputpath_HC_m1);

disp('Loading HC movie2...');
[network_HC_m2, files_HC_m2] = load_network_dir(inputpath_HC_m2);

disp('Loading MCI movie1...');
[network_MCI_m1, files_MCI_m1] = load_network_dir(inputpath_MCI_m1);

disp('Loading MCI movie2...');
[network_MCI_m2, files_MCI_m2] = load_network_dir(inputpath_MCI_m2);

N_HC  = size(network_HC_m1, 1);
N_MCI = size(network_MCI_m1, 1);

[numofchan, ~, numofbands] = size(squeeze(network_HC_m1(1,:,:,:)));

fprintf('HC subjects: %d\n', N_HC);
fprintf('MCI subjects: %d\n', N_MCI);
fprintf('Channels: %d\n', numofchan);
fprintf('Bands: %d\n', numofbands);

assert(size(network_HC_m2,1) == N_HC,  'HC movie1/movie2 subject number mismatch.');
assert(size(network_MCI_m2,1) == N_MCI, 'MCI movie1/movie2 subject number mismatch.');

assert(size(network_HC_m1,2) == numofchan && size(network_HC_m1,3) == numofchan, 'HC movie1 dimension error.');
assert(size(network_HC_m2,2) == numofchan && size(network_HC_m2,3) == numofchan, 'HC movie2 dimension error.');
assert(size(network_MCI_m1,2) == numofchan && size(network_MCI_m1,3) == numofchan, 'MCI movie1 dimension error.');
assert(size(network_MCI_m2,2) == numofchan && size(network_MCI_m2,3) == numofchan, 'MCI movie2 dimension error.');

%% 3. Organize data

D_HC  = network_HC_m2  - network_HC_m1;
D_MCI = network_MCI_m2 - network_MCI_m1;

D_all = cat(1, D_HC, D_MCI);
y_all = [zeros(N_HC,1); ones(N_MCI,1)];

N = N_HC + N_MCI;

D_all(isinf(D_all)) = 0;
D_all(isnan(D_all)) = 0;

%% 4. Prepare output variables

X_train_all = zeros(N, N-1, 10);
y_train_all = zeros(N, N-1);

X_test_all  = zeros(N, 10);
y_test_all  = zeros(N, 1);

train_index_all = zeros(N, N-1);
test_index_all  = zeros(N, 1);

edge_count_negative = zeros(N, length(useBands));
edge_count_positive = zeros(N, length(useBands));
edge_count_candidate = zeros(N, length(useBands));

mask_negative_all = false(N, numofchan, numofchan, length(useBands));
mask_positive_all = false(N, numofchan, numofchan, length(useBands));
mask_candidate_all = false(N, numofchan, numofchan, length(useBands));

upper_mask = triu(true(numofchan), 1);

%% 5. LOOCV nested feature selection

for fold = 1:N

    fprintf('LOOCV fold %d/%d\n', fold, N);

    test_idx = fold;
    train_idx = setdiff(1:N, test_idx);

    test_index_all(fold) = test_idx;
    train_index_all(fold,:) = train_idx;

    y_train = y_all(train_idx);
    y_test  = y_all(test_idx);

    y_train_all(fold,:) = y_train';
    y_test_all(fold) = y_test;

    X_train = zeros(N-1, 10);
    X_test  = zeros(1, 10);

    train_HC_idx  = train_idx(y_all(train_idx) == 0);
    train_MCI_idx = train_idx(y_all(train_idx) == 1);

    for b = 1:length(useBands)

        nb = useBands(b);

        temp_HC  = D_all(train_HC_idx,:,:,nb);
        temp_MCI = D_all(train_MCI_idx,:,:,nb);

        D_HC_vec  = reshape(temp_HC,  length(train_HC_idx),  numofchan*numofchan);
        D_MCI_vec = reshape(temp_MCI, length(train_MCI_idx), numofchan*numofchan);

        [~, p_vec, ~, stats] = ttest2(D_MCI_vec, D_HC_vec);

        t_vec = stats.tstat;

        p_mat = reshape(p_vec, numofchan, numofchan);
        t_mat = reshape(t_vec, numofchan, numofchan);

        p_mat(isnan(p_mat)) = 1;
        t_mat(isnan(t_mat)) = 0;

        candidate_mask = upper_mask & p_mat < alpha;
        sig_mask = make_sig_mask(p_mat, upper_mask, alpha, fdrMode);

        pos_mask = sig_mask & t_mat > 0;
        neg_mask = sig_mask & t_mat < 0;

        mask_candidate_all(fold,:,:,b) = reshape(candidate_mask, [1 numofchan numofchan]);
        mask_negative_all(fold,:,:,b)  = reshape(neg_mask, [1 numofchan numofchan]);
        mask_positive_all(fold,:,:,b)  = reshape(pos_mask, [1 numofchan numofchan]);

        edge_count_candidate(fold,b) = sum(candidate_mask(:));
        edge_count_negative(fold,b)  = sum(neg_mask(:));
        edge_count_positive(fold,b)  = sum(pos_mask(:));

        X_train(:, b)     = extract_network_feature(D_all, train_idx, neg_mask, nb);
        X_train(:, b + 5) = extract_network_feature(D_all, train_idx, pos_mask, nb);

        X_test(1, b)      = extract_network_feature(D_all, test_idx, neg_mask, nb);
        X_test(1, b + 5)  = extract_network_feature(D_all, test_idx, pos_mask, nb);

    end

    X_train_all(fold,:,:) = X_train;
    X_test_all(fold,:) = X_test;

end

%% 6. Save results

save_name = sprintf('LOOCV_nested_features_%s.mat', fdrMode);

save(fullfile(outputpath, save_name), ...
    'X_train_all', 'y_train_all', ...
    'X_test_all', 'y_test_all', ...
    'train_index_all', 'test_index_all', ...
    'edge_count_candidate', ...
    'edge_count_negative', 'edge_count_positive', ...
    'mask_candidate_all', ...
    'mask_negative_all', 'mask_positive_all', ...
    'feature_names', 'band_name', ...
    'N_HC', 'N_MCI', 'N', ...
    'alpha', 'fdrMode', 'useBands', ...
    'files_HC_m1', 'files_HC_m2', ...
    'files_MCI_m1', 'files_MCI_m2', ...
    '-v7.3');

fprintf('Done. Features saved to: %s\n', fullfile(outputpath, save_name));

%% Local functions

function [network_data, files] = load_network_dir(inputpath)
% Load EEG_results.M_zscore_mean from all .mat files in one folder.

    files = dir(fullfile(inputpath, '*.mat'));

    if isempty(files)
        error('No .mat files found in: %s', inputpath);
    end

    [~, idx] = sort({files.name});
    files = files(idx);

    temp = load(fullfile(inputpath, files(1).name));

    if ~isfield(temp, 'EEG_results')
        error('EEG_results not found in file: %s', files(1).name);
    end

    if ~isfield(temp.EEG_results, 'M_zscore_mean')
        error('EEG_results.M_zscore_mean not found in file: %s', files(1).name);
    end

    M = temp.EEG_results.M_zscore_mean;
    [numofchan, ~, numofbands] = size(M);

    Nsub = length(files);
    network_data = zeros(Nsub, numofchan, numofchan, numofbands);

    for i = 1:Nsub
        temp = load(fullfile(inputpath, files(i).name));

        M = temp.EEG_results.M_zscore_mean;
        M(isinf(M)) = 0;
        M(isnan(M)) = 0;

        network_data(i,:,:,:) = M;
    end
end

function sig_mask = make_sig_mask(p_mat, upper_mask, alpha, fdrMode)
% Generate the significant edge mask.

    sig_mask = false(size(p_mat));

    switch lower(fdrMode)

        case 'uncorrected'

            sig_mask = upper_mask & p_mat < alpha;

        case 'candidate_edges_fdr'

            candidate_mask = upper_mask & p_mat < alpha;

            if sum(candidate_mask(:)) > 0

                p_candidate = p_mat(candidate_mask);
                p_candidate(isnan(p_candidate)) = 1;

                q_candidate = mafdr(p_candidate(:), 'BHFDR', true);

                candidate_idx = find(candidate_mask);
                sig_idx = candidate_idx(q_candidate < alpha);

                sig_mask(sig_idx) = true;
            end

        otherwise

            error('Unknown fdrMode: %s', fdrMode);
    end
end

function feat = extract_network_feature(D_all, sub_idx, mask, nb)
% Extract the mean network value inside one mask.

    if isrow(sub_idx)
        sub_idx = sub_idx(:);
    end

    nSub = length(sub_idx);
    feat = zeros(nSub, 1);

    if sum(mask(:)) == 0
        feat(:) = 0;
        return;
    end

    for s = 1:nSub
        mat = squeeze(D_all(sub_idx(s),:,:,nb));
        vals = mat(mask);
        feat(s) = mean(vals, 'omitnan');
    end
end
