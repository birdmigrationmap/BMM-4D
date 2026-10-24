%% Scan cleaning study: select cases for visual review
% Pick a small, diverse review set from the screened volumes:
% - mistnet_more / rho_more: where one method removes much more signal than
%   the other (split by direction so both failure modes are reviewed)
% - both_remove: heavy removal by both methods (do they agree for the right reason?)
% - clean_control: little removal but strong echo, to catch contamination
%   that neither rain filter targets (clutter, insects, wind farms)
% - random: unselected volumes as a sanity check
% Each radar appears at most twice so the set spans different sites.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
data_dir = fullfile(project_root, "prod", "cleaning_study", "data");
load(fullfile(data_dir, "screening.mat"), "screening")
rng(1)

%% Selection controls -----------------------------------------------------
n_mistnet_more = 6;
n_rho_more = 6;
n_both_remove = 4;
n_clean_control = 4;
n_random = 4;
max_per_radar = 2;
both_remove_min = 0.3;
clean_control_max = 0.02;

ok = screening.error == "" & screening.n_gate_reference > 0;
pool = screening(ok, :);
pool.group = strings(height(pool), 1);

%% Select by group ---------------------------------------------------------
% Groups are filled in order; a volume is used only once.
group_names = ["mistnet_more" "rho_more" "both_remove" "clean_control" "random"];
group_sizes = [n_mistnet_more n_rho_more n_both_remove n_clean_control n_random];
for i_g = 1:numel(group_names)
    switch group_names(i_g)
        case "mistnet_more"
            candidate_score = pool.removed_mistnet - pool.removed_rho;
        case "rho_more"
            candidate_score = pool.removed_rho - pool.removed_mistnet;
        case "both_remove"
            candidate_score = min(pool.removed_rho, pool.removed_mistnet);
            candidate_score(candidate_score < both_remove_min) = nan;
        case "clean_control"
            candidate_score = pool.vid_none;
            candidate_score(pool.removed_union > clean_control_max) = nan;
        case "random"
            candidate_score = rand(height(pool), 1);
    end
    candidate_score(pool.group ~= "") = nan;
    [~, order] = sort(candidate_score, "descend", MissingPlacement="last");
    n_selected = 0;
    for i = order(:)'
        if isnan(candidate_score(i)) || n_selected == group_sizes(i_g)
            break
        end
        if nnz(pool.radar(pool.group ~= "") == pool.radar(i)) >= max_per_radar
            continue
        end
        pool.group(i) = group_names(i_g);
        n_selected = n_selected + 1;
    end
end

%% Build the review table -------------------------------------------------
% Review columns are left empty for the visual verdicts.
cases = pool(pool.group ~= "", :);
cases = sortrows(cases, ["group" "disagreement"], ["ascend" "descend"]);
cases.case_id = compose("C%02d", (1:height(cases))');
cases = movevars(cases, ["case_id" "group"], Before=1);
cases.verdict_rho = strings(height(cases), 1);
cases.verdict_mistnet = strings(height(cases), 1);
cases.contamination = strings(height(cases), 1);
cases.notes = strings(height(cases), 1);

save(fullfile(data_dir, "cases.mat"), "cases")
writetable(cases(:, ["case_id" "group" "radar" "scan_time" "season" "removed_rho" "removed_mistnet" ...
    "disagreement" "vid_none" "vid_rho" "vid_mistnet" "verdict_rho" "verdict_mistnet" "contamination" "notes"]), ...
    fullfile(data_dir, "cases_review.csv"))
disp(groupcounts(cases, "group"))
