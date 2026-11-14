function [F, corr] = fobj(cfg, bbar, X)
%FOBJ  式(12)：JPTA 与期望波束的匹配度 Fobj，取值 [0,1]。
%
%   Fobj = sum_k omega_k*|bbar_k' * T*P*d_k| / sum_k omega_k
%
%   由 Cauchy-Schwarz，因 ||bbar_k||=||TPd_k||=1，故 Fobj <= 1。
%   任何 >1 的结果都是 bug（可当单元测试）。

corr = sum(conj(bbar).*X, 1).';              % K x 1，bbar_k' * X_k
F    = sum(cfg.omega.*abs(corr))/sum(cfg.omega);
end
