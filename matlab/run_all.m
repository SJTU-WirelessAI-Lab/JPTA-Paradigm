function run_all(which)
%RUN_ALL  Reproduce the JPTA 2022 results, one entry point.
%
%   run_all            % self-test + all six figures
%   run_all('selftest')
%   run_all('fig4')    % array-gain heat maps
%   run_all('fig5')    % Fobj vs N            (~1.5 min)
%   run_all('fig6')    % Fobj vs kappa        (~1.5 min)
%   run_all('fig7')    % convergence          (~13 min, 100 realizations)
%   run_all('fig9')    % HBF array-gain maps  (~2 min)
%   run_all('fig10')   % Fobj vs N_RF         (~4 min)
%
%   Every figure caches its data in out/*.mat; re-running with an existing
%   cache only redraws. Delete the .mat (or pass true to the fig function)
%   to force a recomputation.
%
%   Console output is ASCII on purpose: MATLAB -batch on a Chinese Windows
%   console prints UTF-8 into a GBK code page and the text is lost.

if nargin < 1 || isempty(which)
    which = 'all';
end

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'core'), fullfile(here,'alg'), fullfile(here,'figs'));

t00 = tic;
switch lower(which)
    case 'selftest'
        selftest();
    case 'fig4'
        fig4_arraygain(cfg_default('N',64,'kappa',64,'maxiter',10));
    case 'fig5'
        fig5_fobj_vs_N(cfg_default());
    case 'fig6'
        fig6_fobj_vs_kappa(cfg_default());
    case 'fig7'
        fig7_convergence(cfg_default(), 100);
    case 'fig9'
        fig9_hbf_arraygain(cfg_default());
    case 'fig10'
        fig10_fobj_vs_nrf(cfg_default());
    case 'all'
        selftest();
        fig4_arraygain(cfg_default('N',64,'kappa',64,'maxiter',10));
        fig5_fobj_vs_N(cfg_default());
        fig6_fobj_vs_kappa(cfg_default());
        fig7_convergence(cfg_default(), 100);
        fig9_hbf_arraygain(cfg_default());
        fig10_fobj_vs_nrf(cfg_default());
    otherwise
        error('run_all: unknown target "%s"', which);
end
fprintf('\n[run_all] total %.1f s\n', toc(t00));
end
