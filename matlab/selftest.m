function selftest()
%SELFTEST  Consistency checks for the forward model and Algorithm 1.
%
%   Checks:
%     1. ||a_k(theta)||^2 = M
%     2. ||T*P*d_k|| = 1
%     3. Fobj of the ideal beam is exactly 1
%     4. Fobj <= 1 for arbitrary parameters (Cauchy-Schwarz)
%     5. The FFT line search of (8a) really sits on the peak
%     6. Algorithm 1 at N=M, kappa=M reproduces the level of Fig. 5(a)/6(a)
%        (behavior 1 ~0.93, behavior 2 ~0.86)
%
%   NOTE: console output is ASCII on purpose -- MATLAB -batch on a Chinese
%   Windows console prints UTF-8 into a GBK code page and the text is lost.
%
%   Usage: cd matlab; selftest

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'core'), fullfile(here,'alg'), fullfile(here,'figs'));
cfg  = cfg_default();

fprintf('\n===== JPTA self-test (M=%d, N=%d, K=%d, kappa=%g) =====\n', ...
        cfg.M, cfg.N, cfg.K, cfg.kappa);

% 1 --------------------------------------------------------------------
A  = array_resp(cfg, 0.3);
e1 = max(abs(sum(abs(A).^2,1) - cfg.M));
report('||||a_k(theta)||^2 - M||_inf', e1, 1e-10);

% 2 --------------------------------------------------------------------
X0 = jpta_forward(cfg, zeros(cfg.N,1), zeros(cfg.M,1));
e2 = max(abs(sum(abs(X0).^2,1) - 1));
report('|||T*P*d_k| - 1||_inf', e2, 1e-10);

% 3 / 4 ----------------------------------------------------------------
[~, bbar1] = gen_beam(cfg, 1);
[~, bbar2] = gen_beam(cfg, 2);
Fid = fobj(cfg, bbar1, bbar1);
report('|Fobj(ideal beam) - 1|', abs(Fid-1), 1e-10);

rng(1);
Fmax = 0;
for t = 1:50
    tau = rand(cfg.N,1)*cfg.tau_max;
    phi = 2*pi*rand(cfg.M,1) - pi;
    Fmax = max(Fmax, fobj(cfg, bbar1, jpta_forward(cfg, tau, phi)));
end
report('max Fobj - 1 (random params, must be <= 0)', max(Fmax-1,0), 0);

% 5 --------------------------------------------------------------------
aa = zeros(cfg.K,1);
[tau_ls, info] = solve_tau_ls(cfg, bbar1, aa);
o0   = direct_obj(cfg, bbar1, aa, tau_ls);
dp   = 80e-12;                                 % 80 ps ~ 13 grid steps
o_up = direct_obj(cfg, bbar1, aa, tau_ls + dp);
o_dn = direct_obj(cfg, bbar1, aa, tau_ls - dp);
worse = sum(o_up(:) > o0(:) + 1e-9) + sum(o_dn(:) > o0(:) + 1e-9);
report('eq.(8a) peak: TTDs improved by 80ps perturbation', worse, 0);
fprintf('      (FFT grid step %.2f ps, with parabolic interpolation)\n', ...
        info.grid_step*1e12);

% 6 --------------------------------------------------------------------
fprintf('\n--- Algorithm 1: N=M=%d, kappa=%g ---\n', cfg.N, cfg.kappa);
for b = 1:2
    [~, bb] = gen_beam(cfg, b);
    t0 = tic;
    [tau, phi, ang, hist] = alg1_jpta(cfg, bb, 'ls');
    el = toc(t0);
    X  = jpta_forward(cfg, tau, phi);
    fprintf('  behavior %d : Fobj = %.4f   time %5.1f s   traj %s\n', ...
            b, fobj(cfg,bb,X), el, mat2str(round(hist.',3)));
end

fprintf('\n--- Algorithm 1 with Lemma 1 (behavior 1 only) ---\n');
[~, bb] = gen_beam(cfg, 1);
t0 = tic; [tau, phi] = alg1_jpta(cfg, bb, 'lemma1'); el = toc(t0);
fprintf('  Lemma 1 : Fobj = %.4f   time %5.1f s\n', ...
        fobj(cfg,bb,jpta_forward(cfg,tau,phi)), el);

fprintf('\n===== self-test done =====\n');
end

% ---------------------------------------------------------------------
function o = direct_obj(cfg, bbar, ang_alpha, tau)
% Direct evaluation of the (8a) objective per TTD (no FFT), used to
% verify that the returned tau really sits on the peak.
c = (cfg.omega.*exp(1j*ang_alpha)).';
B = conj(bbar).*c;                                 % M x K
E = exp(-1j*2*pi*(cfg.f(:).').*tau(cfg.gidx));     % M x K
S = sum(B.*E, 2);                                  % M x 1
o = zeros(cfg.N,1);
for n = 1:cfg.N
    o(n) = sum(abs(S(cfg.gidx==n)));
end
end

function report(name, val, tol)
if val <= tol
    fprintf('  [PASS] %-50s %.3e (tol %.0e)\n', name, val, tol);
else
    fprintf('  [FAIL] %-50s %.3e (tol %.0e)\n', name, val, tol);
end
end
