function [FRF, FBB] = hbf_fc(B, NRF, niter, nrestart, bbar)
%HBF_FC  全连接混合波束成形（式 13/14）的 PE-AltMin 求解。
%
%   [FRF, FBB] = hbf_fc(B, NRF, niter, nrestart, bbar)
%
%   min_{FRF,FBB} || B - FRF*FBB ||_F^2   s.t. |[FRF]_{i,j}| = 1, ||FRF*FBB||_F^2 = Psum
%
%   B    : M x K 期望波束（列 = 各子载波的 b_k，已含功率 sqrt(Psum/(M*K))）
%   NRF  : RF 链数
%   niter: 每次重启的交替迭代次数（默认 200）
%   nrestart : 随机重启次数（默认 8，取残差最小者）
%
%   交替两步（PE-AltMin，见 [11]）：
%     给定 FRF：FBB = pinv(FRF)*B            （最小二乘，无约束）
%     给定 FBB：FRF = exp(j*angle(B*FBB'))   （逐元素取相位，即恒模投影）
%   目标函数单调不增，故固定重启次数即足够。
%
%   注意 B 的列已含功率，||B||_F^2 = Psum，因此收敛后 ||FRF*FBB||_F^2 ≈ Psum
%   自动满足；若需要严格相等，末了对 FBB 再做一次整体缩放。
%
%   ⚠ 选择准则：给了 bbar 时用 Fobj 选优，否则用 Frobenius 残差。
%   NRF 较小时残差准则会掉进局部最优（行为1、NRF=1 时残差最优解只给 Fobj=0.11，
%   比 PC 还低，而 FC 的可行集包含 PC，逻辑上不可能）。原因同样是恒模逼近的
%   局部极小与 Fobj 不单调对应。

if nargin < 3 || isempty(niter),    niter = 200;   end
if nargin < 4 || isempty(nrestart), nrestart = 8;  end
if nargin < 5, bbar = []; end

[M, K] = size(B);
bestScore = -inf; bestRes = inf; FRF = []; FBB = [];

% 确定性初值：B 的前 NRF 个左奇异向量的相位，通常已很好
[Ub,~,~] = svd(B, 'econ');
FRF0 = exp(1j*angle(Ub(:,1:min(NRF,size(Ub,2)))));
if size(FRF0,2) < NRF
    FRF0 = [FRF0, exp(1j*2*pi*rand(M, NRF-size(FRF0,2)))];
end

for t = 1:nrestart
    if t == 1
        F = FRF0;
    else
        F = exp(1j*2*pi*rand(M, NRF));      % 随机相位重启
    end
    for it = 1:niter
        G   = pinv(F)*B;                    % NRF x K
        F   = exp(1j*angle(B*G'));          % 恒模投影
    end
    G   = pinv(F)*B;
    res = norm(B - F*G, 'fro')^2;
    if ~isempty(bbar)
        X = F*G; X = X./vecnorm(X);
        score = sum(abs(sum(conj(bbar).*X, 1)))/size(bbar,2);
        take  = (score > bestScore + 1e-9);
    else
        take  = (res < bestRes);
    end
    if take
        bestScore = score;
        bestRes   = min(bestRes, res);
        FRF = F; FBB = G;
    end
end
end
