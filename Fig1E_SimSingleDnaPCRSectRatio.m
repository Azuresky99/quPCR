clear all; clc
sect = [38.353, 45;
		37.095, 38.07;
		36.523, 37.095;
		36.12, 36.523;
		35.805, 36.12;
		35.547, 35.805;
		35.2, 35.547;
		35.0, 35.2];

E = 1.959 - 1;
fx = 0.12;
threshold = (1+E)^38.34;			% 38.34，make CqP=38.21
MR = fx / threshold;
fprintf("RNAseP MR = %7.6g\n", MR);
nTube = 100000;
nInitDnaCopy = [1, 2, 3, 4, 5, 6, 7, 8];
nStageRound = 20 - round(log(nInitDnaCopy/nInitDnaCopy(1))/log(2));
profileWidth = 1;

for CopyRound = 1:length(nInitDnaCopy)
	fprintf("\n%g. 初始拷贝数 %g\n", CopyRound, nInitDnaCopy(CopyRound));
	nStage1 = nStageRound(CopyRound);
	nDNACopy = ones(nStage1, nTube) * nInitDnaCopy(CopyRound);
	needPlotStage1 = true;
	Ct = zeros(1, nTube);
	
	% 1st stage
	tic
	wbar = waitbar(0,'第一阶段仿真...');			% waitbar
	for roundNum = 2:nStage1
		for tt = 1:nTube
			n = nDNACopy(roundNum-1, tt);
			r = (rand(1, n) <= E);
			p = r + 1;
			n = sum(p);
			nDNACopy(roundNum,tt) = n;
			if mod(tt, 10000)==0
				waitbar(((roundNum-2)*nTube+tt)/((nStage1-1)*nTube+3), wbar, sprintf('第一阶段仿真，第%g轮%g孔，...', roundNum, tt));
			end
		end
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
	
	% 2nd Stage
	wbar = waitbar(1, wbar, '第二阶段仿真...');	% waitbar显示进度条
	for tt = 1:nTube
		Ct(tt) = log(threshold/nDNACopy(nStage1,tt)) / log(E+1) + nStage1 - 1;
	end

	figure, hist(Ct, 1000)
	[counts,centers] = hist(Ct, 1000); title(sprintf('%g拷贝, Ct值分布，总共%d个复孔', nInitDnaCopy(CopyRound), nTube));
	[idx5, idx95] = find90(centers, counts);
	dCenter(CopyRound) = centers(2)-centers(1);
	slideWin = round((centers(idx95)-centers(idx5)) / dCenter(CopyRound) / 10);
	slideWin = floor(slideWin/2) * 2 + 1;
	extCount = slideWin;
	extZero = zeros(1, extCount);
	countsExt = [extZero, counts, extZero];
	centersExtHead = (centers(1)-extCount*dCenter(CopyRound)) : dCenter(CopyRound) : (centers(1)-dCenter(CopyRound));
	centersExtTail = (centers(end)+dCenter(CopyRound)) : dCenter(CopyRound) : (centers(end)+extCount*dCenter(CopyRound));
	centersExt = [centersExtHead, centers, centersExtTail];
	nStart = ceil(extCount/2);
	nEnd = length(countsExt) - extCount;
	ntemplate = normpdf(1:slideWin, nStart, slideWin/5);
	halfWidth = floor(slideWin / 2);
	filtedCounts = zeros(1, length(counts) + extCount + 1);
	for ii = nStart:nEnd
		filtedCounts(ii-halfWidth) = sum(countsExt((ii-halfWidth):(ii+halfWidth)) .* ntemplate);
	end
	filtedCenters = centersExt(nStart:nEnd+nStart);
	figure, bar(centers, counts), hold on, plot(filtedCenters, filtedCounts, 'r', 'linewidth', profileWidth), title([num2str(nInitDnaCopy(CopyRound)), '拷贝，Ct值分布 + 滑窗滤波后的Ct值分布']), legend("Ct", "smooth Ct")
	figure, bar(filtedCenters, filtedCounts), title([num2str(nInitDnaCopy(CopyRound)), '拷贝，滑窗滤波后的Ct值分布'])
	dcmObj = datacursormode;					% Turn on data cursors and return the data cursor mode object
	set(dcmObj, 'UpdateFcn', @myfunction);		% Set the data cursor mode object update function so it uses updateFcn.m
	Ctx = log(threshold/nInitDnaCopy(CopyRound))/log(E+1);
	hold on, plot([Ctx Ctx], [0 max(filtedCounts)*1.2])
	fprintf('\n%g拷贝, Ct理想值 = %g\nCt分布均值 = %g\n以均值代替理想值，误差 %g‰\n', nInitDnaCopy(CopyRound), Ctx, mean(Ct), abs(mean(Ct) - Ctx) / Ctx * 1000);
	PeakIndex = find(filtedCounts == max(filtedCounts));
	fprintf('%g拷贝, 峰值位置 = %g\n\n', nInitDnaCopy(CopyRound), centers(PeakIndex(1)))
	close(wbar);								% waitbar(1, bar, 'sim finished!');
	toc
	RoundLength(CopyRound) = length(filtedCenters);
	centersRound(CopyRound, 1:RoundLength(CopyRound)) = filtedCenters;
	countsRound(CopyRound, :) = counts;
	filtedCountsRound(CopyRound, 1:RoundLength(CopyRound)) = filtedCounts;
	figure(99), hold on, plot(centersRound(CopyRound, 1:RoundLength(CopyRound)), filtedCountsRound(CopyRound, 1:RoundLength(CopyRound)) * dCenter(1)/dCenter(CopyRound));	% dCenter(1)/dCenter(CopyRound)项是各分布高度归一化因子。因为初始模板数越多，结果分布就越集中，每个bin宽度就越小，结果画在同一张图上时，不同初始模板数条件下，衡量y轴高度的统计口径就不一样了。用归一化因子使各种初始模板数条件下的统计口径都一样。
end

CalcFrrAndFar(sect, nInitDnaCopy, centersRound, filtedCountsRound);

function output_txt = myfunction(~, event_obj)
pos = get(event_obj, 'Position');
output_txt = {['X: ', num2str(pos(1), 6)], ...
              ['Y: ', num2str(pos(2), 6)]};
if length(pos) > 2
    output_txt{end+1} = ['Z: ', num2str(pos(3), 6)];
end
end

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