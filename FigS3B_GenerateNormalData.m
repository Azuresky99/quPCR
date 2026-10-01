%% 正态分布的效率值数据生成程序
% NormalData = GenerateNormalData(E, half_width, n)
% 目标：生成 n 个服从该正态分布的数据点
% 输入参数：
%	E,	平均效率。默认值=1.9593
%	half_width,	效率波动边界范围：均值 ± 0.05，即 [1.9093, 2.0093]。该边界对应 95% 置信区间（显著性水平 α = 0.05）
%	n,	生成的正态分布数据的数据量
%	confidence_level, 置信水平，默认99%
% 返回值：	NormalData, 指定的正态分布样本数据。输出返回值时，不自动绘图。
% 示例：
%	NormalData = GenerateNormalData;	% 使用默认参数生成正态分布数组。
%	NormalData = GenerateNormalData(1.9593, 0.05, 100000, 0.99);
%
% 碧云天书，2026年6月24日22:40:02

function NormalData = FigS3B_GenerateNormalData(E, half_width, n, confidence_level)

%% 1. 设置默认参数
if nargin < 1
	E = 1 + 0.9593;				% 均值 (E)
end
if nargin < 2
	half_width = 0.05;			% 边界半宽度
end
if nargin < 3
	n = 100000;					% 指定数据数量
end
if nargin < 4
	confidence_level = 0.99;        % 置信水平 99%
end
lower_bound = E - half_width;   % 下边界 = 1.9093
upper_bound = E + half_width;   % 上边界 = 2.0093
alpha = 1 - confidence_level;   % 显著性水平 1%

if nargout == 0
	fprintf('已知均值 μ = %.4f\n', E);
	fprintf('%g%% 置信区间边界: [%.4f, %.4f]\n', confidence_level*100, lower_bound, upper_bound);
end

%% 2. 根据置信区间反推标准差 sigma
% 标准正态分布的双侧 95% 分位数 z_0.025 ≈ 1.96
% 公式: E ± z_alpha/2 * sigma = 边界
% 因此: sigma = (边界 - E) / z_alpha/2

z = norminv(1 - alpha/2);	% 标准正态分布的双侧 95% 分位数 ≈ 1.96
sigma = half_width / z;		% 计算标准差

if nargout == 0
	fprintf('标准正态分位数 z(0.025) = %.4f\n', z);
	fprintf('反推得到的标准差 σ = %.6f\n', sigma);
end

%% 3. 生成 n 个正态分布随机数
% 使用 randn 生成标准正态分布，再变换为均值为 mu，标准差为 sigma 的正态分布
NormalData = E + sigma * randn(n, 1);

if nargout == 0
	fprintf('\n已生成 %d 个服从 N(%.4f, %.6f²) 的数据点\n', n, E, sigma);
end

if nargout == 0
%% 4. 验证生成数据的统计特性（可选）
	data_mean = mean(NormalData);          % 样本均值
	data_std = std(NormalData);            % 样本标准差
	data_ci = prctile(NormalData, [2.5, 97.5]);  % 样本的 95% 经验分位数区间
	
	fprintf('\n--- 验证统计量 ---\n');
	fprintf('样本均值 = %.6f (理论值 %.4f)\n', data_mean, E);
	fprintf('样本标准差 = %.6f (理论值 %.6f)\n', data_std, sigma);
	fprintf('样本 %g%% 分位数区间: [%.6f, %.6f]\n', confidence_level*100, data_ci(1), data_ci(2));
	fprintf('理论 %g%% 置信区间: [%.4f, %.4f]\n', confidence_level*100, lower_bound, upper_bound);
	
	% 计算实际落在理论区间内的比例
	in_interval = sum(NormalData >= lower_bound & NormalData <= upper_bound) / n * 100;
	fprintf('落在理论区间 [%.4f, %.4f] 内的数据占比: %.2f%%\n', ...
        lower_bound, upper_bound, in_interval);

%% 5. 绘制直方图与理论密度曲线对比（可视化）
	figure;
	histogram(NormalData, 30, 'Normalization', 'pdf', 'FaceColor', [0.7, 0.8, 1]);
	hold on;
	
	% 绘制理论正态密度曲线
	x = linspace(E - 4*sigma, E + 4*sigma, 200);
	y = normpdf(x, E, sigma);
	plot(x, y, 'r-', 'LineWidth', 2);
	
	% 标注均值和置信区间边界
	xline(E, 'k--', '均值', 'LabelOrientation', 'horizontal');
	xline(lower_bound, 'g--', '下边界');
	xline(upper_bound, 'g--', '上边界');
	
	xlabel('数据值');
	ylabel('概率密度');
	title(sprintf('正态分布 N(μ=%.4f, σ=%.6f) 样本直方图 (n=%d)', E, sigma, n));
	legend('样本直方图', '理论密度曲线', '均值', sprintf('%g%% 置信区间边界', confidence_level*100));
	grid on;
	hold off;
end
