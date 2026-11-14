function fig5_fobj_vs_N(cfg0, recompute)
%FIG5_FOBJ_VS_N  图 5：Fobj 随 TTD 数 N 的变化（行为 1 与行为 2）。
%
%   fig5_fobj_vs_N()             % 有缓存就只重绘，没有就重算
%   fig5_fobj_vs_N(cfg, true)    % 强制重算
%
%   四条曲线（照原文）：
%     Heuristic: Algo2 / Algo3      蓝方块实线
%     Algo1 with eqn(8)             红叉实线
%     Algo1 with Lemma1             绿空心圆点线
%     Upper bound                   黑虚线（Cauchy-Schwarz 上界 = 1）
%
%   扫描 N in {1,2,4,8,16,32,64}，kappa 固定为 M（第七节），maxiter=10。
%   N 必须整除 M，否则映射矩阵 P 的等分结构不成立。

if nargin < 1, cfg0 = []; end
if nargin < 2 || isempty(recompute), recompute = false; end
if isempty(cfg0), cfg0 = cfg_default(); end

outdir  = fullfile(fileparts(mfilename('fullpath')),'..','out');
matfile = fullfile(outdir,'fig5_fobj_vs_N.mat');

if exist(matfile,'file') && ~recompute
    S = load(matfile);
    Ns = S.Ns; F = S.F; TAU = S.TAU; %#ok<NASGU>
    fprintf('[fig5] loaded cached data (%s)\n', matfile);
else
    Ns  = [1 2 4 8 16 32 64];
    nN  = numel(Ns);
    F   = nan(nN,3,2);      % [N x {algo1-8a, algo1-lemma1, heuristic} x behavior]
    TAU = cell(1,2);

    t0 = tic;
    for beh = 1:2
        for i = 1:nN
            cfg = cfg_default('N', Ns(i), 'kappa', cfg0.M, 'maxiter', cfg0.maxiter);
            [~, bbar] = gen_beam(cfg, beh);

            % --- Algorithm 1 with (8a) ---
            [tau, phi] = alg1_jpta(cfg, bbar, 'ls');
            F(i,1,beh) = fobj(cfg, bbar, jpta_forward(cfg, tau, phi));

            % --- Algorithm 1 with Lemma 1 ---
            [tau, phi] = alg1_jpta(cfg, bbar, 'lemma1');
            F(i,2,beh) = fobj(cfg, bbar, jpta_forward(cfg, tau, phi));

            % --- heuristic: Algorithm 2 (behavior 1) / Algorithm 3 (behavior 2) ---
            if beh == 1
                [tau, phi] = alg2_heur1(cfg);
            else
                [tau, phi] = alg3_heur2(cfg);
            end
            F(i,3,beh) = fobj(cfg, bbar, jpta_forward(cfg, tau, phi));

            if i == nN, TAU{beh} = tau; end
            fprintf('  beh%d N=%3d : 8a=%.4f  lemma1=%.4f  heur=%.4f  (%.1fs)\n', ...
                    beh, Ns(i), F(i,1,beh), F(i,2,beh), F(i,3,beh), toc(t0));
        end
    end
    save(matfile, 'Ns','F','TAU','-v7');
end

% ------------------------------------------------------------------ 绘图
sty = {struct('c',[0 0 0.85],'m','s','ls','-', 'lw',1.6,'ms',8,'mf','none'), ...
       struct('c',[0.85 0 0],'m','x','ls','-', 'lw',1.6,'ms',9,'mf','none'), ...
       struct('c',[0.25 0.6 0.25],'m','o','ls',':', 'lw',1.8,'ms',9,'mf','none')};
labs = {'Heuristic: Algo2','Algo1 with eqn(8)','Algo1 with Lemma1','Upper bound'};
caps = {'(a) Behavior 1 (\theta_0 = \pi/6, \Delta\theta = \pi/4)', ...
        '(b) Behavior 2 (\theta_1 = -\pi/4, \theta_2 = \pi/6)'};

fig = figure('Visible','off','Position',[100 80 620 860],'Color','w');

for beh = 1:2
    ax = subplot(2,1,beh);
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    if beh == 1, labs{1} = 'Heuristic: Algo2'; else, labs{1} = 'Heuristic: Algo3'; end

    h = gobjects(4,1);
    % 数据列： 1 = eqn(8)，2 = Lemma1，3 = heuristic
    % 样式索引：1 = 蓝方块(heuristic)，2 = 红叉(eqn 8)，3 = 绿圈(Lemma1)
    % 两者编号不同，必须用 col2sty 转换，否则图例标签与实际曲线错配。
    col2sty = [2 3 1];
    for j = [3 2 1]                     % 先蓝、再绿圈，最后红叉压在最上
        k = col2sty(j);
        s = sty{k};
        h(k) = plot(ax, Ns, F(:,j,beh), s.ls, 'Color',s.c, 'Marker',s.m, ...
                    'LineWidth',s.lw, 'MarkerSize',s.ms, 'MarkerFaceColor',s.mf);
    end
    h(4) = plot(ax, [Ns(1) Ns(end)], [1 1], 'k--', 'LineWidth',1.8);

    ylim(ax,[0 1]); xlim(ax,[0 64]);
    xticks(ax, 0:10:60); yticks(ax, 0:0.2:1);
    set(ax,'FontSize',11,'TickDir','out','Layer','top');
    xlabel(ax,'Number of TTDs N','FontSize',12);
    ylabel(ax,'Goodness of fit:  F_{obj}','FontSize',12);

    lg = legend(ax, h, labs, 'Location','southeast','FontSize',11);
    lg.Box = 'on';

    % 面板标题画在坐标框下方（与原图一致）
    text(ax, 0.5, -0.19, caps{beh}, 'Units','normalized', ...
         'HorizontalAlignment','center','FontSize',12);
end

print(fig, fullfile(outdir,'fig5_fobj_vs_N.png'), '-dpng','-r200');
close(fig);
fprintf('[fig5] figure written\n');
end
