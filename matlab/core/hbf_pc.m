function [FRF, FBB] = hbf_pc(B, NRF, niter, nrestart, bbar)
%HBF_PC  部分连接混合波束成形（式 14 + 块稀疏约束）的交替最小化求解。
%
%   [FRF, FBB] = hbf_pc(B, NRF, niter, nrestart, bbar)
%
%   与 hbf_fc 相同的两步交替，但 FRF 只在块对角位置非零：
%       [FRF]_{m,n} = 0,  m ~= M_n,   M_n = {(n-1)M/NRF < m <= n*M/NRF}
%   即每个 RF 链只驱动连续的一段 M/NRF 根天线（图 8 的 PC 结构）。
%
%   ⚠ 选择准则：给了 bbar 时用 Fobj 选优，否则用 Frobenius 残差。
%   原因是 PC 结构下两种解（"各子阵指向同一方向" vs "各子阵分担不同方向"）
%   的残差几乎相同、但 Fobj 差 0.2：行为2、NRF=2 时前者 0.51、后者 0.69。
%   只按残差挑会稳定落到差的那个解，这正是图 10(b) 对不上的原因。

if nargin < 3 || isempty(niter),    niter = 200;   end
if nargin < 4 || isempty(nrestart), nrestart = 8;  end
if nargin < 5, bbar = []; end

[M, K] = size(B);
if mod(M, NRF) ~= 0
    error('hbf_pc: NRF=%d 不能整除 M=%d', NRF, M);
end

% 块划分：mask(m,n)=true 表示第 m 根天线属于第 n 个 RF 链
Pmask = zeros(M, NRF);
for n = 1:NRF
    Pmask(( (n-1)*M/NRF + 1 ):( n*M/NRF ), n) = 1;
end

bestScore = -inf; bestRes = inf;
FRF = []; FBB = [];

% 确定性初值：把 B 的前 NRF 个左奇异向量的相位按块填进 FRF。
[Ub,~,~] = svd(B, 'econ');
Uinit = zeros(M, NRF);
for n = 1:NRF
    rows = ((n-1)*M/NRF + 1):(n*M/NRF);
    Uinit(rows,n) = Ub(rows, min(n,size(Ub,2)));
end
F0 = exp(1j*angle(Uinit)).*Pmask;

for t = 1:nrestart
    if t == 1
        F = F0;
    else
        F = exp(1j*2*pi*rand(M, NRF)).*Pmask;
    end
    for it = 1:niter
        G = pinv(F)*B;                      % NRF x K
        U = B*G';                           % M x K -> 取相位后只保留块内元素
        F = exp(1j*angle(U)).*Pmask;
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
