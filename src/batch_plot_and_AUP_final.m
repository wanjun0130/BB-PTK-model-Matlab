% BB-PTK Model – Batch Post-processing Script
% Adapted from: Yu JY, Rosania GR. Pharm Res. 2010;27(3):457–467.
% This script performs high-throughput post-processing of BB-PTK simulation outputs.
% It assumes that compartmental amounts (mol) have already been computed in upstream simulations.
%% -------------------------------------------------------------------------
% Methodological context:
% (1) Stabilized Nernst–Planck factor
% Implemented in upstream simulation model (not in this script).
% (2) Explicit handling of neutral species (z = 0)
% Implemented in upstream simulation model.
% (3) High-throughput dosimetry extraction
% This script performs batch extraction of pulmonary dosimetry metrics:
%   - Whole-lung burden (mol)
%       includes: aEp + tissue compartments
%   - Tissue-only burden (mol)
%       excludes: aEp
%   - Plasma burden (mol)
%   - AUP (area under burden-time curve)
%   - Tmax
% Additional derived metric:
%   - tissue_to_total_AUP_ratio = AUP_tissue / AUP_total
% Lung burden definitions:
%   Whole-lung:
%   Mlung_total = MaEp + MimEp + McEp + Mint + Msm + MimInt + McEd
%   Tissue-only:
%   Mlung_tissue = MimEp + McEp + Mint + Msm + MimInt + McEd
%% Note:
% This script does NOT perform ODE simulation. It processes exported time-course data from the BB-PTK model.
%% Author: Wanjun Zhang
% Year: 2026
% -------------------------------------------------------------------------

clear; clc;

folderPath = fullfile('output');
dose_all   = 400e-9;   % mol (total deposited dose)
Vp         = 5;        % confirmed in TB/AL
 
% ---- output dirs ----
outLungDir      = fullfile(folderPath, 'fig_timecourse');         % lung only (whole-lung)
outPlasmaDir    = fullfile(folderPath, 'fig_plasma');             % plasma conc + burden
outComboDir     = fullfile(folderPath, 'fig_lung_plasma');        % lung+plasma (mol)
outComboSameDir = fullfile(folderPath, 'fig_lung_plasma_sameY');  % optional: same Y (mol)
 
if ~exist(outLungDir,'dir');      mkdir(outLungDir); end
if ~exist(outPlasmaDir,'dir');    mkdir(outPlasmaDir); end
if ~exist(outComboDir,'dir');     mkdir(outComboDir); end
if ~exist(outComboSameDir,'dir'); mkdir(outComboSameDir); end
 
files = dir(fullfile(folderPath, '*_timecourse.csv'));
if isempty(files)
    error('No *_timecourse.csv found in folder.');
end
 
% ---------- summary columns ----------
AUP_results = cell(numel(files), 17);
%  1  timecourse_file
%  2  AUP_lung_mol_h             (whole-lung incl. aEp)
%  3  AUP_lung_pct_h             (whole-lung incl. aEp)
%  4  Pmax_lung_mol              (whole-lung incl. aEp)
%  5  Tmax_lung_h                (whole-lung incl. aEp)
%  6  n_points
%  7  AUCp_mol_h_per_V
%  8  Cmax_p
%  9  Tmax_p_h
% 10  AUP_plasma_mol_h
% 11  Pmax_plasma_mol
% 12  Tmax_plasma_h
% 13  Vp_used
% 14  AUP_lung_tissue_mol_h      (excluding aEp)
% 15  Pmax_lung_tissue_mol
% 16  Tmax_lung_tissue_h
% 17  tissue_to_total_AUP_ratio
 
for k = 1:numel(files)
    tcFile = fullfile(folderPath, files(k).name);
    T = readtable(tcFile);
 
    % ---- required columns ----
    if ~ismember('time_h', T.Properties.VariableNames)
        error('Missing time_h in file: %s', files(k).name);
    end
    time_h = T.time_h;
 
    % ============================================================
    % Lung burden definitions
    % Priority:
    %   1) Use enhanced columns if present
    %   2) Otherwise reconstruct from compartment columns
    % ============================================================
    if ismember('Mlung_total', T.Properties.VariableNames)
        Mlung_total = T.Mlung_total;
    else
        % fallback for old files
        required = {'MaEp','MimEp','McEp','Mint','Msm','MimInt','McEd'};
        miss = required(~ismember(required, T.Properties.VariableNames));
        if ~isempty(miss)
            error('Missing columns for reconstructing Mlung_total in %s: %s', ...
                files(k).name, strjoin(miss, ', '));
        end
        Mlung_total = T.MaEp + T.MimEp + T.McEp + T.Mint + T.Msm + T.MimInt + T.McEd;
    end
 
    if ismember('Mlung_tissue', T.Properties.VariableNames)
        Mlung_tissue = T.Mlung_tissue;
    else
        % fallback for old files
        required = {'MimEp','McEp','Mint','Msm','MimInt','McEd'};
        miss = required(~ismember(required, T.Properties.VariableNames));
        if ~isempty(miss)
            error('Missing columns for reconstructing Mlung_tissue in %s: %s', ...
                files(k).name, strjoin(miss, ', '));
        end
        Mlung_tissue = T.MimEp + T.McEp + T.Mint + T.Msm + T.MimInt + T.McEd;
    end
 
    % ---- plasma amount (mol) and concentration ----
    if ~ismember('Mp', T.Properties.VariableNames)
        error('Missing Mp in file: %s', files(k).name);
    end
    Mp = T.Mp;            % mol (plasma burden)
    Cp = Mp ./ Vp;        % mol / V
 
    % ============================================================
    % Main AUP / AUC metrics (whole-lung incl. aEp)
    % ============================================================
    AUP_lung_mol_h = trapz(time_h, Mlung_total);
    AUP_lung_pct_h = trapz(time_h, (Mlung_total ./ dose_all) * 100);
 
    [Pmax_lung, idxL] = max(Mlung_total);
    Tmax_lung_h = time_h(idxL);
 
    % ---- tissue-only sensitivity metrics ----
    AUP_lung_tissue_mol_h = trapz(time_h, Mlung_tissue);
    [Pmax_lung_tissue, idxLt] = max(Mlung_tissue);
    Tmax_lung_tissue_h = time_h(idxLt);
 
    if AUP_lung_mol_h > 0
        tissue_to_total_AUP_ratio = AUP_lung_tissue_mol_h / AUP_lung_mol_h;
    else
        tissue_to_total_AUP_ratio = NaN;
    end
 
    % ---- plasma metrics ----
    AUCp_mol_h_per_V = trapz(time_h, Cp);
    [Cmax, idxC] = max(Cp);
    Tmax_p_h = time_h(idxC);
 
    AUP_plasma_mol_h = trapz(time_h, Mp);
    [Pmax_plasma, idxP] = max(Mp);
    Tmax_plasma_h = time_h(idxP);
 
    % ---- store summary ----
    AUP_results{k,1}  = files(k).name;
    AUP_results{k,2}  = AUP_lung_mol_h;
    AUP_results{k,3}  = AUP_lung_pct_h;
    AUP_results{k,4}  = Pmax_lung;
    AUP_results{k,5}  = Tmax_lung_h;
    AUP_results{k,6}  = numel(time_h);
    AUP_results{k,7}  = AUCp_mol_h_per_V;
    AUP_results{k,8}  = Cmax;
    AUP_results{k,9}  = Tmax_p_h;
    AUP_results{k,10} = AUP_plasma_mol_h;
    AUP_results{k,11} = Pmax_plasma;
    AUP_results{k,12} = Tmax_plasma_h;
    AUP_results{k,13} = Vp;
    AUP_results{k,14} = AUP_lung_tissue_mol_h;
    AUP_results{k,15} = Pmax_lung_tissue;
    AUP_results{k,16} = Tmax_lung_tissue_h;
    AUP_results{k,17} = tissue_to_total_AUP_ratio;
 
    % ---- avoid log(0) for plotting ----
    eps_mass_l_total = max(1e-30, 1e-12 * max(Mlung_total));
    Mlung_total_log  = max(Mlung_total, eps_mass_l_total);
 
    eps_mass_l_tissue = max(1e-30, 1e-12 * max(Mlung_tissue));
    Mlung_tissue_log  = max(Mlung_tissue, eps_mass_l_tissue);
 
    eps_mass_p = max(1e-30, 1e-12 * max(Mp));
    Mp_log     = max(Mp, eps_mass_p);
 
    eps_cp = max(1e-30, 1e-12 * max(Cp));
    Cp_log = max(Cp, eps_cp);
 
    % ============================================================
    % (1) Lung-only plot (whole-lung primary + tissue-only dashed)
    % ============================================================
    fig = figure('Visible','off');
    semilogy(time_h, Mlung_total_log, '-', 'LineWidth', 1.8); hold on;
    semilogy(time_h, Mlung_tissue_log, '--', 'LineWidth', 1.5); hold off;
    xlabel('Time (h)');
    ylabel('Lung burden (mol)');
    title([strrep(files(k).name, '_', '\_') '  |  Whole-lung vs tissue-only']);
    grid on;
    legend({'Whole-lung burden (incl. aEp)','Tissue-only burden'}, 'Location','best');
    set(findall(fig,'-property','FontName'),'FontName','Times New Roman');
 
    outPng = fullfile(outLungDir, [files(k).name(1:end-4) '_lung.png']);
    exportgraphics(fig, outPng, 'Resolution', 300);
    close(fig);
 
    % ============================================================
    % (2) Plasma-only plots: (a) Cp and (b) Plasma burden Mp
    % ============================================================
    % (2a) Cp
    fig = figure('Visible','off');
    plot(time_h, Cp, '-', 'LineWidth', 1.5);
    xlabel('Time (h)');
    ylabel('Plasma concentration Cp (mol/V)');
    title([strrep(files(k).name, '_', '\_') '  |  Plasma concentration-time']);
    grid on;
    set(findall(fig,'-property','FontName'),'FontName','Times New Roman');
 
    outPng = fullfile(outPlasmaDir, [files(k).name(1:end-4) '_Cp.png']);
    exportgraphics(fig, outPng, 'Resolution', 300);
    close(fig);
 
    % (2b) Mp
    fig = figure('Visible','off');
    semilogy(time_h, Mp_log, '-', 'LineWidth', 1.5);
    xlabel('Time (h)');
    ylabel('Plasma burden Mp (mol)');
    title([strrep(files(k).name, '_', '\_') '  |  Plasma burden-time']);
    grid on;
    set(findall(fig,'-property','FontName'),'FontName','Times New Roman');
 
    outPng = fullfile(outPlasmaDir, [files(k).name(1:end-4) '_Mp.png']);
    exportgraphics(fig, outPng, 'Resolution', 300);
    close(fig);
 
    % ============================================================
    % (3) Combined plot (dual y-axis): whole-lung vs plasma burdens
    % ============================================================
    fig = figure('Visible','off');
 
    yyaxis left
    semilogy(time_h, Mlung_total_log, '-', 'LineWidth', 1.8);
    ylabel('Whole-lung burden (mol)');
 
    yyaxis right
    semilogy(time_h, Mp_log, '--', 'LineWidth', 1.8);
    ylabel('Plasma burden (mol)');
 
    xlabel('Time (h)');
    title([strrep(files(k).name, '_', '\_') '  |  Whole-lung + Plasma burdens']);
    grid on;
 
    set(findall(fig,'-property','FontName'),'FontName','Times New Roman');
    legend({'Whole-lung burden','Plasma burden'}, 'Location','best');
 
    outPng = fullfile(outComboDir, [files(k).name(1:end-4) '_lung_plasma_mol.png']);
    exportgraphics(fig, outPng, 'Resolution', 300);
    close(fig);
 
    % ============================================================
    % (4) Optional: same Y-axis (whole-lung + tissue-only + plasma)
    % ============================================================
    fig = figure('Visible','off');
    semilogy(time_h, Mlung_total_log, '-', 'LineWidth', 1.8); hold on;
    semilogy(time_h, Mlung_tissue_log, '--', 'LineWidth', 1.5);
    semilogy(time_h, Mp_log, ':', 'LineWidth', 1.8); hold off;
    xlabel('Time (h)');
    ylabel('Burden (mol)');
    title([strrep(files(k).name, '_', '\_') '  |  Whole-lung vs tissue-only vs plasma']);
    grid on;
    legend({'Whole-lung burden','Tissue-only burden','Plasma burden'}, 'Location','best');
    set(findall(fig,'-property','FontName'),'FontName','Times New Roman');
 
    outPng = fullfile(outComboSameDir, [files(k).name(1:end-4) '_lung_vs_plasma_sameY.png']);
    exportgraphics(fig, outPng, 'Resolution', 300);
    close(fig);
end
 
% ---- write AUP summary (xlsx + csv) ----
AUP_table = cell2table(AUP_results, 'VariableNames', ...
    {'timecourse_file', ...
     'AUP_lung_mol_h','AUP_lung_pct_h','Pmax_lung_mol','Tmax_lung_h','n_points', ...
     'AUCp_mol_h_per_V','Cmax_p','Tmax_p_h', ...
     'AUP_plasma_mol_h','Pmax_plasma_mol','Tmax_plasma_h', ...
     'Vp_used', ...
     'AUP_lung_tissue_mol_h','Pmax_lung_tissue_mol','Tmax_lung_tissue_h', ...
     'tissue_to_total_AUP_ratio'});
 
outXlsx = fullfile(folderPath, 'AUP_sum.xlsx');
outCsv  = fullfile(folderPath, ' AUP_sum.csv');
 
writetable(AUP_table, outXlsx);
writetable(AUP_table, outCsv);
 
fprintf('Done.\n');
fprintf('Lung-only figs (whole + tissue): %s\n', outLungDir);
fprintf('Plasma figs (Cp + Mp):           %s\n', outPlasmaDir);
fprintf('Combo dual-axis:                 %s\n', outComboDir);
fprintf('Combo sameY:                     %s\n', outComboSameDir);
fprintf('AUP summary XLSX:                %s\n', outXlsx);
fprintf('AUP summary CSV :                %s\n', outCsv);
