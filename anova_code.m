%% Mixed ANOVA for EEG network data
clear; clc;

%% Parameters
N1 = 33;                     % MCI group
N2 = 39;                     % HC group
numofchan = 63;
numofbands = 6;
alpha = 0.05;

band_name = {'delta'; 'theta'; 'alpha'; 'beta'; 'gamma'; 'fullband'};

readOrgDataFlag = true;
inputDataFile = 'orgNetData.mat';
outputDir = 'results';

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Load data
if readOrgDataFlag
    [network1_1, network1_2, network2_1, network2_2, ...
        y1_1, y1_2, y2_1, y2_2, ~] = readANOVAorgData(N1, N2, numofchan, numofbands);
else
    load(inputDataFile, ...
        'network1_1', 'network1_2', 'network2_1', 'network2_2', ...
        'y1_1', 'y1_2', 'y2_1', 'y2_2');
end

%% Mixed ANOVA
X2 = ones(2 * N1 + 2 * N2, 1);
X2(2 * N2 + 1 : 2 * N2 + 2 * N1) = 2;

X3 = ones(2 * N2 + 2 * N1, 1);
X3(N2 + 1 : 2 * N2) = 2;
X3(2 * N2 + N1 + 1 : 2 * N2 + 2 * N1) = 2;

X4 = [1:N2, 1:N2, N2 + 1:N1 + N2, N2 + 1:N1 + N2]';

p_movie = nan(numofchan^2, numofbands);
f_movie = nan(numofchan^2, numofbands);
p_patient = nan(numofchan^2, numofbands);
f_patient = nan(numofchan^2, numofbands);
p_inter = nan(numofchan^2, numofbands);
f_inter = nan(numofchan^2, numofbands);

for nb = 1:numofbands
    for side = 1:numofchan^2
        MCI_movie1 = y1_1(:, side, nb);
        MCI_movie2 = y1_2(:, side, nb);
        HC_movie1 = y2_1(:, side, nb);
        HC_movie2 = y2_2(:, side, nb);

        X = [];
        X(:, 1) = [HC_movie1; HC_movie2; MCI_movie1; MCI_movie2];
        X(:, 2) = X2;
        X(:, 3) = X3;
        X(:, 4) = X4;

        [~, ~, ~, Fs, Ps] = mixed_between_within_anova(X, 1);

        p_movie(side, nb) = Ps{3};
        f_movie(side, nb) = Fs{3};

        p_patient(side, nb) = Ps{1};
        f_patient(side, nb) = Fs{1};

        p_inter(side, nb) = Ps{4};
        f_inter(side, nb) = Fs{4};
    end
end

p_movie = reshape(p_movie, numofchan, numofchan, numofbands);
f_movie = reshape(f_movie, numofchan, numofchan, numofbands);
p_patient = reshape(p_patient, numofchan, numofchan, numofbands);
f_patient = reshape(f_patient, numofchan, numofchan, numofbands);
p_inter = reshape(p_inter, numofchan, numofchan, numofbands);
f_inter = reshape(f_inter, numofchan, numofchan, numofbands);

%% Select effect
flagOfConcern = 'inter';     % 'patient', 'movie', or 'inter'

if strcmp(flagOfConcern, 'patient')
    nameOfConcern = 'patientEffect';
    p_concern = p_patient;
    f_concern = f_patient;
elseif strcmp(flagOfConcern, 'movie')
    nameOfConcern = 'movieEffect';
    p_concern = p_movie;
    f_concern = f_movie;
elseif strcmp(flagOfConcern, 'inter')
    nameOfConcern = 'interEffect';
    p_concern = p_inter;
    f_concern = f_inter;
else
    error('Invalid flagOfConcern.');
end

resultDir = fullfile(outputDir, nameOfConcern);
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

%% Save thresholded F values
Result = struct();

for nb = 1:numofbands - 1
    p_temp = p_concern(:, :, nb);
    f_temp = f_concern(:, :, nb);

    p_temp(isnan(p_temp)) = 0;
    f_temp(isnan(f_temp)) = 0;
    f_temp(p_temp > alpha) = 0;
    f_temp(p_temp == 0) = 0;

    Result.f(:, :, nb) = f_temp;
end

%% Mark significant edges
p_marked = zeros(numofchan, numofchan, numofbands);

for nb = 1:numofbands
    p_temp = p_concern(:, :, nb);

    p_temp(isnan(p_temp)) = 0;
    p_temp(p_temp > alpha) = 0;

    p_marked(:, :, nb) = p_temp;
end

Result.p_marked = p_marked;

%% Select positive or negative interaction edges
flagOfPorN = 1;              % 1 = positive, 0 = negative

if flagOfPorN == 1
    nameOfPorN = 'positive';
else
    nameOfPorN = 'negative';
end

maxEdges = numofchan * (numofchan + 1) / 2;
index_PorN = zeros(maxEdges, numofbands);
num_sidePorN = zeros(1, numofbands);

markedSide1 = struct();
markedSide2 = struct();
line_mci_m1_PorN = struct();
line_mci_m2_PorN = struct();
line_scd_m1_PorN = struct();
line_scd_m2_PorN = struct();

for nb = 1:numofbands - 1
    mci_m1 = network1_1(:, :, :, nb);
    mci_m2 = network1_2(:, :, :, nb);
    scd_m1 = network2_1(:, :, :, nb);
    scd_m2 = network2_2(:, :, :, nb);

    nodeOfMarkedSide1 = [];
    nodeOfMarkedSide2 = [];
    num_PorN = 0;

    mci_m1_PorN = [];
    mci_m2_PorN = [];
    scd_m1_PorN = [];
    scd_m2_PorN = [];

    for b1 = 1:numofchan
        for b2 = b1:numofchan
            if p_marked(b1, b2, nb) ~= 0 && p_marked(b1, b2, nb) <= alpha
                nodeOfMarkedSide1 = [nodeOfMarkedSide1, b1];
                nodeOfMarkedSide2 = [nodeOfMarkedSide2, b2];
            end
        end
    end

    markedSide1.(band_name{nb}) = nodeOfMarkedSide1;
    markedSide2.(band_name{nb}) = nodeOfMarkedSide2;

    for i = 1:length(nodeOfMarkedSide1)
        b1 = nodeOfMarkedSide1(i);
        b2 = nodeOfMarkedSide2(i);

        mci_m1_side = squeeze(mci_m1(:, b1, b2));
        mci_m2_side = squeeze(mci_m2(:, b1, b2));
        scd_m1_side = squeeze(scd_m1(:, b1, b2));
        scd_m2_side = squeeze(scd_m2(:, b1, b2));

        yy1_side = [mean(mci_m1_side), mean(scd_m1_side)];
        yy2_side = [mean(mci_m2_side), mean(scd_m2_side)];

        interactionValue = (yy2_side(1) - yy1_side(1)) - (yy2_side(2) - yy1_side(2));

        if (flagOfPorN == 1 && interactionValue > 0) || ...
                (flagOfPorN ~= 1 && interactionValue < 0)

            num_PorN = num_PorN + 1;
            index_PorN(num_PorN, nb) = i;

            mci_m1_PorN(:, num_PorN) = mci_m1_side;
            mci_m2_PorN(:, num_PorN) = mci_m2_side;
            scd_m1_PorN(:, num_PorN) = scd_m1_side;
            scd_m2_PorN(:, num_PorN) = scd_m2_side;
        end
    end

    num_sidePorN(nb) = num_PorN;

    line_mci_m1_PorN.(band_name{nb}) = mci_m1_PorN;
    line_mci_m2_PorN.(band_name{nb}) = mci_m2_PorN;
    line_scd_m1_PorN.(band_name{nb}) = scd_m1_PorN;
    line_scd_m2_PorN.(band_name{nb}) = scd_m2_PorN;
end

Result.index_PorN = index_PorN;
Result.num_sidePorN = num_sidePorN;
Result.markedSide1 = markedSide1;
Result.markedSide2 = markedSide2;

save(fullfile(resultDir, ['line_', nameOfPorN, '.mat']), ...
    'line_mci_m1_PorN', 'line_mci_m2_PorN', ...
    'line_scd_m1_PorN', 'line_scd_m2_PorN');

%% Post-hoc t-tests on significant ANOVA edges
if strcmp(flagOfConcern, 'patient')
    for nb = 1:numofbands
        MCI_m1 = y1_1(:, :, nb);
        MCI_m2 = y1_2(:, :, nb);
        HC_m1 = y2_1(:, :, nb);
        HC_m2 = y2_2(:, :, nb);

        M = MCI_m1 + MCI_m2;
        N = HC_m1 + HC_m2;

        [~, p, ~, stats] = ttest2(M, N);
        post_t_patient(:, nb) = stats.tstat;
        post_p_patient(:, nb) = p;
    end

    t = reshape(post_t_patient, numofchan, numofchan, numofbands);
    p = reshape(post_p_patient, numofchan, numofchan, numofbands);

elseif strcmp(flagOfConcern, 'movie')
    for nb = 1:numofbands
        MCI_m1 = y1_1(:, :, nb);
        MCI_m2 = y1_2(:, :, nb);
        HC_m1 = y2_1(:, :, nb);
        HC_m2 = y2_2(:, :, nb);

        M = [MCI_m2; HC_m2];
        N = [MCI_m1; HC_m1];

        [~, p, ~, stats] = ttest(M, N);
        post_t_movie(:, nb) = stats.tstat;
        post_p_movie(:, nb) = p;
    end

    t = reshape(post_t_movie, numofchan, numofchan, numofbands);
    p = reshape(post_p_movie, numofchan, numofchan, numofbands);

elseif strcmp(flagOfConcern, 'inter')
    for nb = 1:numofbands
        MCI_m1 = y1_1(:, :, nb);
        MCI_m2 = y1_2(:, :, nb);
        HC_m1 = y2_1(:, :, nb);
        HC_m2 = y2_2(:, :, nb);

        M = MCI_m2 - MCI_m1;
        N = HC_m2 - HC_m1;

        [~, p, ~, stats] = ttest2(M, N);
        post_t_inter(:, nb) = stats.tstat;
        post_p_inter(:, nb) = p;
    end

    t = reshape(post_t_inter, numofchan, numofchan, numofbands);
    p = reshape(post_p_inter, numofchan, numofchan, numofbands);
end

t(p_marked == 0) = 0;
p(p_marked == 0) = 0;

%% FDR correction
fdrFlag = true;

for nb = 1:numofbands - 1
    if fdrFlag
        p_temp = triu(reshape(p(:, :, nb), numofchan, numofchan));
        fdr_temp = p_temp;

        p_for_fdr = p_temp(~isnan(p_temp));
        p_for_fdr = p_for_fdr(p_for_fdr ~= 0)';

        if ~isempty(p_for_fdr)
            fdr = mafdr(p_for_fdr, 'BHFDR', true);
            fdr_temp(~isnan(fdr_temp) & fdr_temp ~= 0) = fdr;
            fdr_temp = fdr_temp + fdr_temp';
            p1(:, :, nb) = fdr_temp;
            p_ = p1(:, :, nb);
        else
            p_ = ones(numofchan);
        end
    else
        p_ = p(:, :, nb);
    end

    p_(isnan(p_)) = 1;

    t_ = t(:, :, nb);
    t_(isnan(t_)) = 0;
    t_(p_ > alpha) = 0;

    Result.p_ttest(:, :, nb) = p_;
    Result.t_ttest(:, :, nb) = t_;

    t_marked = zeros(numofchan);

    for i = 1:size(index_PorN, 1)
        if index_PorN(i, nb) ~= 0
            edgeIndex = index_PorN(i, nb);
            b1_temp = markedSide1.(band_name{nb})(edgeIndex);
            b2_temp = markedSide2.(band_name{nb})(edgeIndex);
            t_marked(b1_temp, b2_temp) = t_(b1_temp, b2_temp);
        end
    end

    Result.t_marked(:, :, nb) = t_marked;
end

save(fullfile(resultDir, ['anova_Results_', nameOfPorN, '.mat']), 'Result');
