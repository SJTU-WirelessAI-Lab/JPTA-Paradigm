function G = array_gain(cfg, X, theta_deg)
%ARRAY_GAIN  图 4 的阵列增益：|a_k(theta)' * X_k|^2，角度 x 子载波，线性标度。
%
%   X 为 M x K 的模拟波束矩阵：传 jpta_forward 的输出得到 JPTA 的增益，
%   传 bbar 得到理想波束的增益。
%
%   因 ||a_k||=sqrt(M)、||X_k||=1，峰值恰为 M（本例 64）——色标是线性 0..M，
%   不是 dB。误当 dB 会以为差了好几个数量级。

th = theta_deg(:)*pi/180;
G  = zeros(numel(th), cfg.K);
m  = (0:cfg.M-1).';
fr = cfg.f(:).'/cfg.f0;

for i = 1:numel(th)
    v = exp(1j*pi*m*(sin(th(i))*fr));        % M x K，a_k(theta)
    G(i,:) = abs(sum(conj(v).*X, 1)).^2;     % 期望峰值 M
end
end
