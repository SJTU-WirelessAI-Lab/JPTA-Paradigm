function fig10_fobj_vs_nrf(cfg, recompute)
%FIG10_FOBJ_VS_NRF  图 10：HBF 的 Fobj 随 RF 链数 N_RF 的变化。
%
%   fig10_fobj_vs_nrf()             % 有缓存就只重绘
%   fig10_fobj_vs_nrf(cfg, true)    % 强制重算
%
%   (a) 行为 1，两个场景：
%         theta0 = pi/6, Dtheta = pi/4;   theta0 = 0, Dtheta = 2pi/3
%       每个场景画 HBF-FC、HBF-PC，另加一条 JPTA（1 RF 链、N=64 TTD）参考线。
%   (b) 行为 2，theta1 = -pi/4, theta2 = pi/6：HBF-FC、HBF-PC、JPTA 参考。
%
%   N_RF 取 {1,2,4,8,16,32,64}（PC 要求整除 M=64）。
%   JPTA 参考与 N_RF 无关，画成水平线。

if nargin < 1 || isempty(cfg), cfg = cfg_default(); end
if nargin < 2 || isempty(recompute), recompute = false; end

outdir  = fullfile(fileparts(mfilename('fullpath')),'..','out');
matfile = fullfile(outdir,'fig10_fobj_vs_nrf.mat');

NRF  = [1 2 4 8 16 32 64];
NRST = 4;      % 重启次数（SVD 初值已很好，4 次够）
NIT  = 150;

% 场景定义
scen = struct( ...
    'tag', {'b1a','b1b','b2'}, ...
    'beh', {1, 1, 2}, ...
    'cfgargs', {{'theta0',pi/6,'dtheta',pi/4}, ...
                {'theta0',0,   'dtheta',2*pi/3}, ...
                {'theta1',-pi/4,'theta2',pi/6}});

if exist(matfile,'file') && ~recompute
    S = load(matfile);
    R = S.R; NRF = S.NRF;
    fprintf('[fig10] loaded cached data (%s)\n', matfile);
else
    R = struct();
    t0 = tic;
    for s = 1:numel(scen)
        c = cfg_default(scen(s).cfgargs{:}, 'N',64, 'kappa',64, 'maxiter',10);
        b = gen_beam(c, scen(s).beh);
        [~, bbar] = gen_beam(c, scen(s).beh);

        ff = nan(numel(NRF),1); fp = nan(numel(NRF),1);
        for i = 1:numel(NRF)
            [F1,G1] = hbf_fc(b, NRF(i), NIT, NRST, bbar);
            X1 = F1*G1; X1 = X1./vecnorm(X1);
            ff(i) = fobj(c, bbar, X1);

            [F2,G2] = hbf_pc(b, NRF(i), NIT, NRST, bbar);
            X2 = F2*G2; X2 = X2./vecnorm(X2);
            fp(i) = fobj(c, bbar, X2);
        end

        % JPTA 参考：1 RF 链 + N=64 TTD
        [tau, phi] = alg1_jpta(c, bbar, 'ls');
        fj = fobj(c, bbar, jpta_forward(c, tau, phi));

        R.(scen(s).tag) = struct('fc',ff,'pc',fp,'jpta',fj);
        fprintf('  %-3s : FC=[%s]  PC=[%s]  JPTA=%.4f  (%.0fs)\n', scen(s).tag, ...
                num2str(ff.', '%.3f '), num2str(fp.', '%.3f '), fj, toc(t0));
    end
    save(matfile, 'R','NRF','-v7');
end

% ------------------------------------------------------------------ 绘图
fig = figure('Visible','off','Position',[100 70 640 880],'Color','w');

% ---------------- (a) Behavior 1 ----------------
ax = subplot(2,1,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
h = gobjects(6,1);
h(1) = plot(ax, NRF, R.b1a.fc, 'r-',  'LineWidth',2.0);
h(2) = plot(ax, NRF, R.b1a.pc, 'b-o', 'LineWidth',1.8,'MarkerSize',8);
h(3) = plot(ax, NRF, R.b1b.fc, 'r--', 'LineWidth',2.0);
h(4) = plot(ax, NRF, R.b1b.pc, 'b--o','LineWidth',1.8,'MarkerSize',8);
h(5) = plot(ax, NRF, R.b1a.jpta*ones(size(NRF)), 'k-x','LineWidth',1.8,'MarkerSize',9);
h(6) = plot(ax, NRF, R.b1b.jpta*ones(size(NRF)), 'k--x','LineWidth',1.8,'MarkerSize',9);
ylim(ax,[0 1]); xlim(ax,[0 64]); xticks(ax, 0:10:60); yticks(ax, 0:0.2:1);
set(ax,'FontSize',11,'TickDir','out','Layer','top');
xlabel(ax,'Number of RF chains N_{RF}','FontSize',12);
ylabel(ax,'Goodness of fit:  F_{obj}','FontSize',12);
lg = legend(ax, h, {'HBF-FC, \theta_0 = \pi/6, \Delta\theta = \pi/4', ...
                    'HBF-PC, \theta_0 = \pi/6, \Delta\theta = \pi/4', ...
                    'HBF-FC, \theta_0 = 0, \Delta\theta = 2\pi/3', ...
                    'HBF-PC, \theta_0 = 0, \Delta\theta = 2\pi/3', ...
                    'JPTA, \theta_0 = \pi/6, \Delta\theta = \pi/4', ...
                    'JPTA, \theta_0 = 0, \Delta\theta = 2\pi/3'}, ...
                   'Location','southeast','FontSize',9);   % 右下空白，避免遮住 NRF>32 的曲线
lg.Box = 'on';
text(ax, 0.5, -0.22, '(a) Behavior 1', 'Units','normalized', ...
     'HorizontalAlignment','center','FontSize',12);

% ---------------- (b) Behavior 2 ----------------
ax = subplot(2,1,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
h = gobjects(3,1);
h(1) = plot(ax, NRF, R.b2.fc, 'r-',  'LineWidth',2.0);
h(2) = plot(ax, NRF, R.b2.pc, 'b-o', 'LineWidth',1.8,'MarkerSize',8);
h(3) = plot(ax, NRF, R.b2.jpta*ones(size(NRF)), 'k-x','LineWidth',1.8,'MarkerSize',9);
% 下界取 0.45 而不是原文的 0.6：NRF=1 时 FC 与 PC 都退化为纯模拟波束，
% 理论上 Fobj≈0.5（单波束只能照顾两个角度之一），若照原文取 0.6 会把该点裁掉。
ylim(ax,[0.45 1.0]); xlim(ax,[0 64]); xticks(ax, 0:10:60); yticks(ax, 0.5:0.1:1.0);
set(ax,'FontSize',11,'TickDir','out','Layer','top');
xlabel(ax,'Number of RF chains N_{RF}','FontSize',12);
ylabel(ax,'Goodness of fit:  F_{obj}','FontSize',12);
lg = legend(ax, h, {'HBF-FC, \theta_1 = -\pi/4, \theta_2 = \pi/6', ...
                    'HBF-PC, \theta_1 = -\pi/4, \theta_2 = \pi/6', ...
                    'JPTA, \theta_1 = -\pi/4, \theta_2 = \pi/6'}, ...
                   'Location','southeast','FontSize',10);
lg.Box = 'on';
text(ax, 0.5, -0.22, '(b) Behavior 2', 'Units','normalized', ...
     'HorizontalAlignment','center','FontSize',12);

print(fig, fullfile(outdir,'fig10_fobj_vs_nrf.png'), '-dpng','-r200');
close(fig);
fprintf('[fig10] figure written\n');
end
