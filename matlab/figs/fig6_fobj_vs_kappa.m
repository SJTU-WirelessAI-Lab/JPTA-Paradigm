function fig6_fobj_vs_kappa(cfg0, recompute)
%FIG6_FOBJ_VS_KAPPA  图 6：Fobj 随 TTD 最大时延参数 kappa 的变化。
%
%   fig6_fobj_vs_kappa()             % 有缓存就只重绘
%   fig6_fobj_vs_kappa(cfg, true)    % 强制重算
%
%   曲线与图 5 相同（4 条），N 固定为 M（第七节），maxiter=10。
%   横轴为对数刻度，kappa 采样点按原图取：
%     (a) 行为 1: kappa in {1,2,4,8,16,32,64}
%     (b) 行为 2: kappa in {0.1,0.2,0.4,0.8,1.6,3.2,6.4}
%   注意 kappa 的单位是"1/W 归一化后的时延跨度"，kappa=M 即 [0,6.4] ns。

if nargin < 1, cfg0 = []; end
if nargin < 2 || isempty(recompute), recompute = false; end
if isempty(cfg0), cfg0 = cfg_default(); end

outdir  = fullfile(fileparts(mfilename('fullpath')),'..','out');
matfile = fullfile(outdir,'fig6_fobj_vs_kappa.mat');

kap  = {[1 2 4 8 16 32 64], [0.1 0.2 0.4 0.8 1.6 3.2 6.4]};
cfgN = [1 2];                          % 每个行为对应的 panel

if exist(matfile,'file') && ~recompute
    S = load(matfile);
    kap = S.kap; F = S.F;
    fprintf('[fig6] loaded cached data (%s)\n', matfile);
else
    F = cell(1,2);
    t0 = tic;
    for p = 1:2
        beh = cfgN(p);
        kv  = kap{p};
        F{p} = nan(numel(kv),3);
        for i = 1:numel(kv)
            cfg = cfg_default('N', cfg0.M, 'kappa', kv(i), 'maxiter', cfg0.maxiter);
            [~, bbar] = gen_beam(cfg, beh);

            [tau, phi] = alg1_jpta(cfg, bbar, 'ls');
            F{p}(i,1) = fobj(cfg, bbar, jpta_forward(cfg, tau, phi));

            [tau, phi] = alg1_jpta(cfg, bbar, 'lemma1');
            F{p}(i,2) = fobj(cfg, bbar, jpta_forward(cfg, tau, phi));

            if beh == 1
                [tau, phi] = alg2_heur1(cfg);
            else
                [tau, phi] = alg3_heur2(cfg);
            end
            F{p}(i,3) = fobj(cfg, bbar, jpta_forward(cfg, tau, phi));

            fprintf('  beh%d kappa=%6.2f : 8a=%.4f  lemma1=%.4f  heur=%.4f  (%.1fs)\n', ...
                    beh, kv(i), F{p}(i,1), F{p}(i,2), F{p}(i,3), toc(t0));
        end
    end
    save(matfile, 'kap','F','-v7');
end

% ------------------------------------------------------------------ 绘图
sty = {struct('c',[0 0 0.85],'m','s','ls','-', 'lw',1.6,'ms',8,'mf','none'), ...
       struct('c',[0.85 0 0],'m','x','ls','-', 'lw',1.6,'ms',9,'mf','none'), ...
       struct('c',[0.25 0.6 0.25],'m','o','ls',':', 'lw',1.8,'ms',9,'mf','none')};
caps = {'(a) Behavior 1 (\theta_0 = \pi/6, \Delta\theta = \pi/4)', ...
        '(b) Behavior 2 (\theta_1 = -\pi/4, \theta_2 = \pi/6)'};

fig = figure('Visible','off','Position',[100 80 620 860],'Color','w');

for p = 1:2
    ax = subplot(2,1,p);
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    kv = kap{p};

    labs = {'Heuristic: Algo2','Algo1 with eqn(8)','Algo1 with Lemma1','Upper bound'};
    if p == 2, labs{1} = 'Heuristic: Algo3'; end

    h = gobjects(4,1);
    % 数据列： 1 = eqn(8)，2 = Lemma1，3 = heuristic
    % 样式索引：1 = 蓝方块(heuristic)，2 = 红叉(eqn 8)，3 = 绿圈(Lemma1)
    col2sty = [2 3 1];
    for j = [3 2 1]                     % 先蓝、再绿圈，最后红叉压在最上
        k = col2sty(j);
        s = sty{k};
        h(k) = semilogx(ax, kv, F{p}(:,j), s.ls, 'Color',s.c, 'Marker',s.m, ...
                        'LineWidth',s.lw, 'MarkerSize',s.ms, 'MarkerFaceColor',s.mf);
    end

    % x 范围照原图留一个 decade 的余量（(a) 10^0..10^2，(b) 10^-1..10^1），
    % 否则最右端的采样点会贴着边框、且 10^2 的刻度不显示。
    if p == 1
        xr = [1 100];   xt = [1 10 100];
    else
        xr = [0.1 10];  xt = [0.1 1 10];
    end
    xlim(ax, xr); xticks(ax, xt);
    ylim(ax, [0 1]);
    yticks(ax, 0:0.2:1);
    h(4) = plot(ax, xr, [1 1], 'k--', 'LineWidth',1.8);
    set(ax,'FontSize',11,'TickDir','out','Layer','top','XScale','log','XMinorTick','on');
    xlabel(ax,'TTD max. delay parameter:  \kappa','FontSize',12);
    ylabel(ax,'Goodness of fit:  F_{obj}','FontSize',12);

    lg = legend(ax, h, labs, 'Location','southeast','FontSize',11);
    lg.Box = 'on';

    text(ax, 0.5, -0.19, caps{p}, 'Units','normalized', ...
         'HorizontalAlignment','center','FontSize',12);
end

print(fig, fullfile(outdir,'fig6_fobj_vs_kappa.png'), '-dpng','-r200');
close(fig);
fprintf('[fig6] figure written\n');
end
