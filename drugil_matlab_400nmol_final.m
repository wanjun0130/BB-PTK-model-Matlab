% Adapted from: Yu JY, Rosania GR. Cell-based multiscale computational modeling of small molecule absorption and retention in the lungs. Pharm Res. 2010;27(3):457–467.
% This implementation is an independent adaptation and extension of the published framework, rather than the original code released by the authors.
% Major methodological extensions and improvements:
% (1) Stabilized Nernst–Planck factor
% Original form: phi(N) = N / (exp(N) - 1)
% Modified with a numerically stable implementation:
% phi(N) ≈ 1 - N/2 + N^2/12 for |N| < 1e-6
% phi(N) = N / (exp(N) - 1) otherwise
% This avoids numerical instability when N is close to zero.
% (2) Explicit handling of neutral species (z = 0)
% Neutral molecules bypass ionization partitioning in the implementation.
% The neutral-only branch sets ionized fractions to zero and avoids non-physical ionization-related transport artifacts.
% (3) Extension to a high-throughput BB-PTK platform
% This implementation extends the original transport framework into a burden-based pulmonary tissue kinetics (BB-PTK) platform by:
% converting concentration-based outputs to compartmental amounts using compartment volumes
% enabling automated batch processing across multiple compounds
% supporting automated extraction of dosimetry-oriented outputs such as lung burden (mol), plasma burden (mol), Tmax, and AUP
% Lung burden definition: Z_lung = Z(1)+Z(2)+Z(3)+Z(4)+Z(5)+Z(6)+Z(7), where Z denotes compartmental amount after converting concentration-based state variables using compartment volumes.
% Author: Wanjun Zhang
% Year: 2026
% -------------------------------------------------------------------------
% Robust full version (enhanced):
% 1) Align TB/AL time grids by interpolation onto a common t_common (union)
% 2) Pad matrices to same number of compartments before adding
% 3) Correct T50 logic: plasma mass (compartment 8) reaches 50% of total dose
% 4) Export per-compound full time-course CSV
% 5) Export per-compound summary CSV (header + 1 row)
% 6) Export summary_all.csv for all compounds
% 7) Export merged Excel file
% 8) Whole-lung burden used for model¨Cdata comparison includes aEp
% 9) Summary additionally reports Mlung_total / Mlung_tissue / Clung_total / Clung_tissue
 
clear; clc;
 
%% ---------------- USER SETTINGS ----------------
inputXlsx = fullfile('input', 'input_template.xlsx');
folderPath = fullfile('output'); 
dose_all = 400e-9;
dose_tb  = 280e-9;
dose_al  = dose_all - dose_tb;
 
use_last_if_no_T50 = true;
 
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
 
%% ---------------- READ INPUT ----------------
dataTable = readtable(inputXlsx, 'PreserveVariableNames', true);
 
% ---- Robust numeric extraction ----
logPN_ls = col_to_double(dataTable, 'logPN');
pKa_ls   = col_to_double(dataTable, 'pKa');
z_ls     = col_to_double(dataTable, 'z');
 
% ---- sanitize z to numeric -1/0/+1 ----
z_ls(~isfinite(z_ls)) = 0;
z_ls(abs(z_ls) < 1e-9) = 0;
z_ls(z_ls > 0) = 1;
z_ls(z_ls < 0) = -1;
 
DrugName_ls = dataTable{:, 'DrugName'};
obs_t_ls    = dataTable{:, 'obs_T50'};
 
%% ---------------- MERGE VAR NAMES ----------------
mergeVarNames = {'DrugName','logP','pKa','obs_T50','T50_min','t50_tb_min','t50_al_min', ...
    'aEp_pct','lung_pct','int_pct','plasma_pct', ...
    'Mlung_total','Mlung_tissue','Clung_total','Clung_tissue', ...
    'CaEp','Clung','Cint','Cp', ...
    'CaEpN','CaEpD','CaEpT','CintN','CintD','CintT', ...
    'PLung_h_tb','PSL_h_tb','PLung_h_al','PSL_h_al', ...
    'PMlung_max_tb','TLung_max_tb','PMlung_max_al','TLung_max_al', ...
    'PMimEp_max_tb','PMcEp_max_tb','PMint_max_tb','PMsm_max_tb','PMimInt_max_tb','PMcEd_max_tb', ...
    'PMimEp_max_al','PMcEp_max_al','PMint_max_al','PMsm_max_al','PMimInt_max_al','PMcEd_max_al'};
 
numVars = numel(mergeVarNames);
results = {};
 
%% ---------------- MASTER SUMMARY FILE ----------------
summaryAllPath = fullfile(folderPath, 'summary_all.csv');
fidAll = fopen(summaryAllPath, 'w');
if fidAll == -1
    error('Cannot open summary_all.csv for writing.');
end
fprintf(fidAll, '%s\n', strjoin(mergeVarNames, ','));
 
%% ---------------- MAIN LOOP ----------------
for i = 1:length(logPN_ls)
 
    DrugName = char(DrugName_ls(i));
    logPN    = logPN_ls(i);
    pKa      = pKa_ls(i);
    z        = z_ls(i);
    obs_t    = obs_t_ls(i);
 
    fprintf('Processing %d/%d: %s\n', i, length(logPN_ls), DrugName);
 
    try
        % ---------- HARD INPUT CHECKS ----------
        if ~isfinite(logPN)
            error('Input logPN is NaN/Inf (cannot run).');
        end
 
        if z ~= 0 && ~isfinite(pKa)
            error('Input pKa is NaN/Inf for charged molecule (z=%g).', z);
        end
 
        % ---- Run TB ----
        [t_half_a_tb, ~, t_tb, Y_tb, M_v_tb, M_v_n_tb, M_v_d_tb, M_v_t_tb, ...
            PLung_h_tb, PSL_h_tb, PMlung_max_tb, TLung_max_tb, ...
            PMimEp_max_tb, PMcEp_max_tb, PMint_max_tb, PMsm_max_tb, PMimInt_max_tb, PMcEd_max_tb] = ...
            LungModel_Rat_tb_D_15_dense_huge_free_final(logPN, pKa, z, dose_tb);
 
        % ---- Run AL ----
        [t_half_a_al, ~, t_al, Y_al, M_v_al, M_v_n_al, M_v_d_al, M_v_al_t, ...
            PLung_h_al, PSL_h_al, PMlung_max_al, TLung_max_al, ...
            PMimEp_max_al, PMcEp_max_al, PMint_max_al, PMsm_max_al, PMimInt_max_al, PMcEd_max_al] = ...
            LungModel_Rat_al_D_15_dense_huge_free_final(logPN, pKa, z, dose_al);
 
        if isempty(Y_tb) || isempty(Y_al) || isempty(t_tb) || isempty(t_al)
            fprintf('Empty output for %s, skipping.\n', DrugName);
            continue;
        end
 
        % ---- Common time grid ----
        t_common = unique([t_tb; t_al]);
 
        Y_tb_i = interp1(t_tb, Y_tb, t_common, 'pchip', 'extrap');
        Y_al_i = interp1(t_al, Y_al, t_common, 'pchip', 'extrap');
        Y_tb_i = max(Y_tb_i, 0);
        Y_al_i = max(Y_al_i, 0);
 
        % ---- Convert to mass ----
        Z_tb   = Y_tb_i * M_v_tb;
        Z_n_tb = Y_tb_i * M_v_n_tb;
        Z_d_tb = Y_tb_i * M_v_d_tb;
        Z_t_tb = Y_tb_i * M_v_t_tb;
 
        Z_al   = Y_al_i * M_v_al;
        Z_n_al = Y_al_i * M_v_n_al;
        Z_d_al = Y_al_i * M_v_d_al;
        Z_t_al = Y_al_i * M_v_al_t;
 
        % ---- Pad to same number of compartments ----
        ncol = max(size(Z_tb,2), size(Z_al,2));
        if ncol < 8
            fprintf('Unexpected compartment number (%d) for %s, skipping.\n', ncol, DrugName);
            continue;
        end
 
        Z_tb   = pad_cols(Z_tb,   ncol);
        Z_n_tb = pad_cols(Z_n_tb, ncol);
        Z_d_tb = pad_cols(Z_d_tb, ncol);
        Z_t_tb = pad_cols(Z_t_tb, ncol);
 
        Z_al   = pad_cols(Z_al,   ncol);
        Z_n_al = pad_cols(Z_n_al, ncol);
        Z_d_al = pad_cols(Z_d_al, ncol);
        Z_t_al = pad_cols(Z_t_al, ncol);
 
        % ---- Combine TB + AL ----
        Z_all   = Z_tb   + Z_al;
        Z_n_all = Z_n_tb + Z_n_al;
        Z_d_all = Z_d_tb + Z_d_al;
        Z_t_all = Z_t_tb + Z_t_al;
 
        % ---- Volumes ----
        M_v_tb_p = pad_square(M_v_tb, ncol);
        M_v_al_p = pad_square(M_v_al, ncol);
        V = diag(M_v_tb_p + M_v_al_p)';   % 1 x ncol
 
        % ---- Indices ----
        idx_aEp    = 1;
        idx_imEp   = 2;
        idx_cEp    = 3;
        idx_int    = 4;
        idx_sm     = 5;
        idx_imInt  = 6;
        idx_cEd    = 7;
        idx_plasma = 8;
 
        % ---- T50: plasma reaches 50% of total dose ----
        j = find(Z_all(:, idx_plasma) >= 0.5 * dose_all, 1, 'first');
        if isempty(j)
            if use_last_if_no_T50
                j = length(t_common);
                fprintf('T50 not reached for %s; using last time point.\n', DrugName);
            else
                fprintf('T50 not reached for %s; skipping.\n', DrugName);
                continue;
            end
        end
 
        t50_min = t_common(j) / 60;
 
        % ---- Mass at time j ----
        MaEp   = Z_all(j, idx_aEp);
        MimEp  = Z_all(j, idx_imEp);
        McEp   = Z_all(j, idx_cEp);
        Mint   = Z_all(j, idx_int);
        Msm    = Z_all(j, idx_sm);
        MimInt = Z_all(j, idx_imInt);
        McEd   = Z_all(j, idx_cEd);
        Mp     = Z_all(j, idx_plasma);
 
        % ---- Lung burden definitions ----
        Mlung_tissue = MimEp + McEp + Mint + Msm + MimInt + McEd;
        Mlung_total  = MaEp + MimEp + McEp + Mint + Msm + MimInt + McEd;
 
        % ---- Volumes for concentration calculation ----
        Vlung_tissue = V(idx_imEp) + V(idx_cEp) + V(idx_int) + V(idx_sm) + V(idx_imInt) + V(idx_cEd);
        Vlung_total  = V(idx_aEp) + V(idx_imEp) + V(idx_cEp) + V(idx_int) + V(idx_sm) + V(idx_imInt) + V(idx_cEd);
 
        % ---- Concentrations ----
        CaEp = MaEp / V(idx_aEp);
        Cint = Mint / V(idx_int);
        Cp   = Mp   / V(idx_plasma);
 
        Clung_total  = Mlung_total  / Vlung_total;
        Clung_tissue = Mlung_tissue / Vlung_tissue;
 
        % backward-compatible alias
        Clung = Clung_total;
 
        % ---- Percent dose ----
        aEp_pct    = MaEp        / dose_all * 100;
        lung_pct   = Mlung_total / dose_all * 100;
        int_pct    = Mint        / dose_all * 100;
        plasma_pct = Mp          / dose_all * 100;
 
        % ---- Free concentrations ----
        CaEpN = Z_n_all(j, idx_aEp) / V(idx_aEp);
        CaEpD = Z_d_all(j, idx_aEp) / V(idx_aEp);
        CaEpT = Z_t_all(j, idx_aEp) / V(idx_aEp);
 
        CintN = Z_n_all(j, idx_int) / V(idx_int);
        CintD = Z_d_all(j, idx_int) / V(idx_int);
        CintT = Z_t_all(j, idx_int) / V(idx_int);
 
        % ---- Export full time-course CSV ----
        tcFile = fullfile(folderPath, sprintf('%d_timecourse.csv', i));
        Mlung_total_ts  = sum(Z_all(:,1:7), 2);
        Mlung_tissue_ts = sum(Z_all(:,2:7), 2);
 
        Vlung_total_ts  = V(idx_aEp) + V(idx_imEp) + V(idx_cEp) + V(idx_int) + V(idx_sm) + V(idx_imInt) + V(idx_cEd);
        Vlung_tissue_ts = V(idx_imEp) + V(idx_cEp) + V(idx_int) + V(idx_sm) + V(idx_imInt) + V(idx_cEd);
 
        Clung_total_ts  = Mlung_total_ts  ./ Vlung_total_ts;
        Clung_tissue_ts = Mlung_tissue_ts ./ Vlung_tissue_ts;
 
        tcTable = table( ...
            t_common/3600, ...
            Z_all(:,1), Z_all(:,2), Z_all(:,3), Z_all(:,4), Z_all(:,5), Z_all(:,6), Z_all(:,7), Z_all(:,8), ...
            Mlung_total_ts, Mlung_tissue_ts, Clung_total_ts, Clung_tissue_ts, ...
            'VariableNames', {'time_h','MaEp','MimEp','McEp','Mint','Msm','MimInt','McEd','Mp', ...
                              'Mlung_total','Mlung_tissue','Clung_total','Clung_tissue'} ...
        );
        writetable(tcTable, tcFile);
 
        % ---- Build FIXED-LENGTH results row ----
        row = cell(1, numVars);
        row(:) = {NaN};
 
        row{1}  = DrugName;
        row{2}  = logPN;
        row{3}  = pKa;
        row{4}  = obs_t;
        row{5}  = t50_min;
        row{6}  = t_half_a_tb/60;
        row{7}  = t_half_a_al/60;
 
        row{8}  = aEp_pct;
        row{9}  = lung_pct;
        row{10} = int_pct;
        row{11} = plasma_pct;
 
        row{12} = Mlung_total;
        row{13} = Mlung_tissue;
        row{14} = Clung_total;
        row{15} = Clung_tissue;
 
        row{16} = CaEp;
        row{17} = Clung;
        row{18} = Cint;
        row{19} = Cp;
 
        row{20} = CaEpN;
        row{21} = CaEpD;
        row{22} = CaEpT;
        row{23} = CintN;
        row{24} = CintD;
        row{25} = CintT;
 
        row{26} = PLung_h_tb;
        row{27} = PSL_h_tb;
        row{28} = PLung_h_al;
        row{29} = PSL_h_al;
 
        row{30} = PMlung_max_tb;
        row{31} = TLung_max_tb;
        row{32} = PMlung_max_al;
        row{33} = TLung_max_al;
 
        row{34} = PMimEp_max_tb;
        row{35} = PMcEp_max_tb;
        row{36} = PMint_max_tb;
        row{37} = PMsm_max_tb;
        row{38} = PMimInt_max_tb;
        row{39} = PMcEd_max_tb;
 
        row{40} = PMimEp_max_al;
        row{41} = PMcEp_max_al;
        row{42} = PMint_max_al;
        row{43} = PMsm_max_al;
        row{44} = PMimInt_max_al;
        row{45} = PMcEd_max_al;
 
        % ---- Save per-compound summary CSV ----
        perSummaryPath = fullfile(folderPath, sprintf('%d_summary.csv', i));
        fidOne = fopen(perSummaryPath, 'w');
        if fidOne ~= -1
            fprintf(fidOne, '%s\n', strjoin(mergeVarNames, ','));
            write_row_csv(fidOne, row);
            fclose(fidOne);
        else
            fprintf('Warning: cannot write per-compound summary for %s\n', DrugName);
        end
 
        % ---- Append to master summary_all.csv ----
        write_row_csv(fidAll, row);
 
        % ---- Keep in memory for Excel ----
        results = [results; row];
 
    catch ME
        fprintf('Error processing %s (i=%d): %s\n', DrugName, i, ME.message);
        for ss = 1:min(numel(ME.stack),3)
            fprintf('  at %s (line %d)\n', ME.stack(ss).name, ME.stack(ss).line);
        end
        fprintf('  DEBUG inputs: logPN=%s, pKa=%s, z=%s\n', mat2str(logPN), mat2str(pKa), mat2str(z));
        continue;
    end
end
 
%% ---------------- CLOSE MASTER SUMMARY ----------------
fclose(fidAll);
 
%% ---------------- MERGE RESULTS TO EXCEL ----------------
if isempty(results)
    error('No successful compounds were processed. Please check model functions and inputs.');
end
 
mergeFilePath = fullfile(folderPath, 'summary_merged.xlsx');
mergeTable = cell2table(results, 'VariableNames', mergeVarNames);
writetable(mergeTable, mergeFilePath);
 
fprintf('Done.\n');
fprintf('1) summary_all.csv: %s\n', summaryAllPath);
fprintf('2) merged excel:    %s\n', mergeFilePath);
fprintf('3) per-compound timecourse + summary in folder:\n   %s\n', folderPath);
 
%% ---------------- HELPER FUNCTIONS ----------------
function A2 = pad_cols(A, ncol)
    if size(A,2) < ncol
        A2 = [A, zeros(size(A,1), ncol - size(A,2))];
    else
        A2 = A(:,1:ncol);
    end
end
 
function M2 = pad_square(M, n)
    [r,c] = size(M);
    if r==n && c==n
        M2 = M;
        return;
    end
    M2 = zeros(n,n);
    rr = min(r,n);
    cc = min(c,n);
    M2(1:rr,1:cc) = M(1:rr,1:cc);
end
 
function write_row_csv(fid, row)
    N = numel(row);
    for c = 1:N
        val = row{c};
        if isstring(val); val = char(val); end
 
        if ischar(val)
            val = strrep(val, '"', '""');
            fprintf(fid, '"%s"', val);
        elseif isnumeric(val)
            if isempty(val) || numel(val)~=1 || isnan(val)
                fprintf(fid, '');
            else
                fprintf(fid, '%.10g', val);
            end
        else
            fprintf(fid, '');
        end
 
        if c < N
            fprintf(fid, ',');
        else
            fprintf(fid, '\n');
        end
    end
end
 
function x = col_to_double(T, varName)
    raw = T{:, varName};
 
    if iscell(raw)
        x = nan(size(raw));
        for ii = 1:numel(raw)
            v = raw{ii};
            if isempty(v)
                x(ii) = NaN;
            elseif isstring(v) || ischar(v)
                x(ii) = str2double(strtrim(char(v)));
            else
                x(ii) = double(v);
            end
        end
    elseif iscategorical(raw)
        x = str2double(strtrim(string(raw)));
    elseif isstring(raw)
        x = str2double(strtrim(raw));
    else
        x = double(raw);
    end
end
