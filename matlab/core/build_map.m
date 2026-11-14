function P = build_map(M, N)
%BUILD_MAP  式(1) 的映射矩阵 P（M x N），每行恰好一个非零元。
%
%   论文第七节：M_n = {m : (n-1)M/N < m <= nM/N}，即连续分块连接，
%   相邻天线尽量共用同一个 TTD，便于布线（图 3）。

m = (1:M).';
P = zeros(M, N);
for n = 1:N
    P(:,n) = ((n-1)*M/N < m) & (m <= n*M/N);
end

if any(sum(P,2) ~= 1)
    error('build_map: 映射矩阵每行必须恰好一个非零元');
end
end
