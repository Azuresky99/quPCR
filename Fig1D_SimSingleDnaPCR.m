% simulation of replication
% the program is base on Matlab

clear all, clc
E = 1.959 - 1;
threshold = (1+E)^38.34;
nWell = 1000000;	
nInitDnaCopy = 1;
nStageRound = 20;

for CopyRound = 1:length(nInitDnaCopy)
	fprintf("\n%g. 初始拷贝数 %g\n", CopyRound, nInitDnaCopy(CopyRound));
	nStage1 = nStageRound(CopyRound);
	nDNACopy = ones(nStage1, nWell) * nInitDnaCopy(CopyRound);
	needPlotStage1 = true;
	Cq = zeros(1, nWell);
	
	% 1st stage
	tic
	wbar = waitbar(0,'第一阶段仿真...');	
	for roundNum = 2:nStage1
		prevRow = nDNACopy(roundNum-1, :);
		parfor tt = 1:nWell	
			n = prevRow(tt);
			r = (rand(1, n) <= E);
			p = r + 1;
			n = sum(p);
			nDNACopy(roundNum,tt) = n;
		end
		waitbar(((roundNum-2)*nWell)/((nStage1-1)*nWell), wbar, sprintf('第一阶段仿真，第%g轮，...', roundNum));
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
	wbar = waitbar(1, wbar, '第二阶段仿真...');		% waitbar show progress
	parfor tt = 1:nWell
		Cq(tt) = log(threshold/nDNACopy(nStage1,tt)) / log(E+1) + nStage1 - 1;
	end

	% show results
	figure, hist(Cq, 1000)
	[counts,centers] = hist(Cq, 1000); title(sprintf('%g拷贝, Cq值分布，总共%d个复孔', nInitDnaCopy(CopyRound), nWell));
	[idx5, idx95] = find90(centers, counts);
	slideWin = round(0.15 / (centers(2) - centers(1)));
	slideWin = floor(slideWin/2) * 2 + 1;
	nStart = ceil(slideWin / 2);
	nEnd = length(centers) - floor(slideWin / 2);
	slideCounts = counts;
	halfWidth = floor(slideWin / 2);
	ntemplate = normpdf(1:slideWin, nStart, slideWin/5);
	for ii = nStart:nEnd
		slideCounts(ii) = sum(counts((ii-halfWidth):(ii+halfWidth)) .* ntemplate);
	end
	figure, bar(centers, counts), hold on, plot(centers, slideCounts, 'r'), title([num2str(nInitDnaCopy(CopyRound)), '拷贝，Cq值分布 + 滑窗滤波后的Cq值分布']), legend("Cq", "smooth Cq")
	figure, bar(centers, slideCounts), title([num2str(nInitDnaCopy(CopyRound)), '拷贝，滑窗滤波后的Cq值分布'])
	dcmObj = datacursormode;  % Turn on data cursors and return the
	set(dcmObj, 'UpdateFcn', @myfunction);  % Set the data cursor mode object update
	Cqx = log(threshold/nInitDnaCopy(CopyRound))/log(E+1);
	hold on, plot([Cqx Cqx], [0 max(slideCounts)*1.2])
	fprintf('\n%g拷贝, Cq理想值 = %g\nCq分布均值 = %g\n以均值代替理想值，误差 %g‰\n', nInitDnaCopy(CopyRound), Cqx, mean(Cq), abs(mean(Cq) - Cqx) / Cqx * 1000);
	PeakIndex = find(slideCounts(nStart:nEnd) == max(slideCounts(nStart:nEnd)));
	fprintf('%g拷贝, 峰值位置 = %g\n\n', nInitDnaCopy(CopyRound), centers(PeakIndex(1)+nStart-1))
	close(wbar);
	toc
end

function output_txt = myfunction(obj, event_obj)
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