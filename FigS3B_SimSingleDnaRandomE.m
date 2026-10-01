E = 1.9593 - 1;					% 平均效率
devE = 0.2;						% 允许的效率波动范围
confidence_level = 0.95;		% 置信水平95%
fx = 0.12;						% Cq阈值
threshold = (1+E)^38.335;		% 手工凑阈值，使Cq峰值=38.21。阈值就是拷贝数，单位：条
MR = fx / threshold;
fprintf("RNAseP MR = %7.6g\n", MR);
nTube = 1000000;				% 反应复孔的数量。100万个耗时很长
nInitDnaCopy = 1;				% 初始状态下，复孔内的DNA拷贝数
nStageRound = 20;				% 对应不同的nInitDnaCpoy初始拷贝数，第一阶段完成时的轮次

for CopyRound = 1:length(nInitDnaCopy)
	fprintf("\n%g. 初始拷贝数 %g\n", CopyRound, nInitDnaCopy(CopyRound));
	nStage1 = nStageRound(CopyRound);
	nDNACopy = ones(nStage1, nTube) * nInitDnaCopy(CopyRound);
	NormalE = reshape(FigS3B_GenerateNormalData(E+1, devE/2, nStage1*nTube, confidence_level), [nStage1 nTube]);
	needPlotStage1 = true;
	Cq = zeros(1, nTube);
	
	wbar = waitbar(0,'第一阶段仿真...');
	for roundNum = 2:nStage1			
    	prevRow = nDNACopy(roundNum - 1, :);
    	curRow = zeros(1, nTube);			
    	parfor tt = 1:nTube
        	n = prevRow(tt);
			r = (rand(1, n) <= (NormalE(roundNum, tt)-1));
			p = r + 1;
			n = sum(p);
        	curRow(tt) = n;
    	end
    	nDNACopy(roundNum, :) = curRow;
    	waitbar((roundNum - 1) / (nStage1 - 1), wbar, sprintf('第一阶段仿真，第 %g 轮完成', roundNum));
	end
	if needPlotStage1
		waitbar(1, wbar, '计算概率密度分布图...');
		ymin = 1;									
		ymax = max(nDNACopy(end, :));
		xi = 1:0.01:nStage1;
		yi = interp1(1:nStage1, single(nDNACopy), xi, 'spline');
		hmap = zeros(length(xi), 500);				
		hmap100 = zeros(length(xi), 500);			
		ybins = ymin: (ymax-1) / (500-1) : ymax;
		for ii = 1:length(xi)
			hh = hist(yi(ii, :), ybins);
			hmap(ii,:) = log10(hh);
			hmap100(ii,:) = hmap(ii,:) / max(hmap(ii,:));
		end
		figure, imagesc([1 nStage1], [1 ymax], hmap'), axis xy, cbar=colorbar; colormap('jet'), xlabel('轮次'), ylabel('拷贝数'), title('扩增反应的概率密度分布图'); cbar.Title.String = '复孔数量';	
		figure, imagesc([1 nStage1], [1 ymax], hmap100'), axis xy, cbar=colorbar; colormap('jet'), xlabel('轮次'), ylabel('拷贝数'), title('扩增反应的归一化概率密度分布图'); cbar.Title.String = '归一化的复孔占比';
	end
	
	wbar = waitbar(1, wbar, '第二阶段仿真...');
	for tt = 1:nTube
		Cq(tt) = log(threshold/nDNACopy(nStage1,tt)) / log(E+1) + nStage1 - 1;
	end

	figure, hist(Cq, 1000)
	[counts,centers] = hist(Cq, 1000); title(sprintf('%g拷贝, Cq值分布，总共%d个复孔', nInitDnaCopy(CopyRound), nTube));
	[idx5, idx95] = find90(centers, counts);
	slideWin = round((centers(idx95)-centers(idx5)) / 10 / (centers(2) - centers(1)));
	slideWin = floor(slideWin/2) * 2 + 1;			% 确保滑窗大小为奇数
	nStart = ceil(slideWin / 2);
	nEnd = length(centers) - floor(slideWin / 2);
	slideCounts = counts;
	halfWidth = floor(slideWin / 2);
	ntemplate = normpdf(1:slideWin, nStart, slideWin/5);
	for ii = nStart:nEnd
		slideCounts(ii) = sum(counts((ii-halfWidth):(ii+halfWidth)) .* ntemplate);	% 高斯滑窗
	end
	figure, bar(centers, counts, 'BarWidth', 1), hold on, plot(centers(nStart:nEnd), slideCounts(nStart:nEnd), 'r'), title(sprintf('%d拷贝, E = %g±%g, 置信水平 = %g%%, Cq值分布 + 滑窗滤波后的Cq值分布', nInitDnaCopy(CopyRound), E+1, devE/2, confidence_level*100)), legend("Cq", "smooth Cq");
	Cqx = log(threshold/nInitDnaCopy(CopyRound))/log(E+1);
	hold on, plot([Cqx Cqx], [0 max(slideCounts)*1.2])
	figure, bar(centers, slideCounts, 'BarWidth', 1), title(sprintf('%d拷贝, E = %g±%g, 置信水平 = %g%%, 滑窗滤波后的Cq值分布', nInitDnaCopy(CopyRound), E+1, devE/2, confidence_level*100));
	dcmObj = datacursormode;
	set(dcmObj, 'UpdateFcn', @myfunction);
	Cqx = log(threshold/nInitDnaCopy(CopyRound))/log(E+1);
	hold on, plot([Cqx Cqx], [0 max(slideCounts)*1.2])
	fprintf('\n%g拷贝, Cq理想值 = %g\nCq分布均值 = %g\n以均值代替理想值，误差 %g‰\n', nInitDnaCopy(CopyRound), Cqx, mean(Cq), abs(mean(Cq) - Cqx) / Cqx * 1000);
	PeakIndex = find(slideCounts(nStart:nEnd) == max(slideCounts(nStart:nEnd)));
	fprintf('%g拷贝, 峰值位置 = %g\n\n', nInitDnaCopy(CopyRound), centers(PeakIndex(1)+nStart-1))
	close(wbar);
end

function output_txt = myfunction(obj, event_obj)
	pos = get(event_obj, 'Position');
	output_txt = {['X: ', num2str(pos(1), 6)], ...
              	['Y: ', num2str(pos(2), 6)]};
	if length(pos) > 2
    	output_txt{end+1} = ['Z: ', num2str(pos(3), 6)];
	end
end

% 找分布曲线90%数据区域
function [idx5, idx95] = find90(x, y)
	area = trapz(x, y);
	if abs(area - 1) > 0.01
		y = y / area;
	end
	cum_prob = cumsum(y);
	cum_prob = cum_prob / max(cum_prob);
	idx5 = find(cum_prob >= 0.05, 1);
	idx95 = find(cum_prob >= 0.95, 1);
end