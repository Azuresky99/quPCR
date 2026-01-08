% the simulation is divide into two stage. In the first stage, every DNA is
% replicated according the probability (AE) in each cycle, until the total
% DNA number in a well acheive the copy number which correspond to that
% of 1 DNA replicate 20 cycles. In the second stage, the Cq value of every
% well is directly calculated according the expenatial amplification
% equation based on the DNA copy number when the first stage finished.
% 
% 2025-12-28 18:37:18

clear all; clc
E = 1.9596 - 1;					% Amplify efficiency
fx = 0.12;						% Cq threshold
threshold = (1+E)^38.335;		% set the copy number while amplicons cross the Cq, and make the Peak of Cq profile equal to 38.21
MR = fx / threshold;			% fluorescence efficiency
nTube = 1000000;				% the number of replication wells for each initial copy number
nInitDnaCopy = [1, 2, 4, 8, 16, 100, 1000, 10000];	% the initial copy number of wells
nStageRound = [20, 19, 18, 17, 16, 14, 11, 7];		% the end cycles of the first stage for different initial copy number

for CopyRound = 1:length(nInitDnaCopy)
	fprintf("\n%g. Initial copy number: %g\n", CopyRound, nInitDnaCopy(CopyRound));
	nStage1 = nStageRound(CopyRound);
	nDNACopy = ones(nStage1, nTube) * nInitDnaCopy(CopyRound);
	needPlotStage1 = true;
	Cq = zeros(1, nTube);
	
	% The first stage
	tic
	wbar = waitbar(0,'First stage...');
	for roundNum = 2:nStage1
		for tt = 1:nTube
			n = nDNACopy(roundNum-1, tt);
			r = (rand(1, n) <= E);
			p = r + 1;
			n = sum(p);
			nDNACopy(roundNum,tt) = n;
			if mod(tt, 10000)==0
				waitbar(((roundNum-2)*nTube+tt)/((nStage1-1)*nTube+3), wbar, sprintf('First stage, the %g round %g well，...', roundNum, tt));
			end
		end
	end
	
	% The second stage
	wbar = waitbar(1, wbar, 'Second stage...');		% progress bar
	for tt = 1:nTube
		Cq(tt) = log(threshold/nDNACopy(nStage1,tt)) / log(E+1) + nStage1 - 1;
	end

	% plot
	figure, hist(Cq, 1000)
	[counts,centers] = hist(Cq, 1000); title(sprintf('%g copy, Cq value spread，%d wells in total', nInitDnaCopy(CopyRound), nTube));
	% denoise by slide window
	[idx5, idx95] = find90(centers, counts);
	slideWin = round((centers(idx95)-centers(idx5)) / 10 / (centers(2) - centers(1)));
	slideWin = floor(slideWin/2) * 2 + 1;
	nStart = ceil(slideWin / 2);
	nEnd = length(centers) - floor(slideWin / 2);
	slideCounts = counts;
	halfWidth = floor(slideWin / 2);
	ntemplate = normpdf(1:slideWin, nStart, slideWin/5);					% Gaussian curve
	for ii = nStart:nEnd
		slideCounts(ii) = sum(counts((ii-halfWidth):(ii+halfWidth)) .* ntemplate);
	end
	figure, bar(centers, counts), hold on, plot(centers, slideCounts, 'r'), title([num2str(nInitDnaCopy(CopyRound)), 'Copy, Cq value spread + Cq value after sliding filter']), legend("Cq", "smooth Cq")
	figure, bar(centers, slideCounts), title([num2str(nInitDnaCopy(CopyRound)), 'Copy, Cq value spread after sliding filter'])
	% set high precision cursor
	dcmObj = datacursormode;% Turn on data cursors and return the
							%   data cursor mode object
	set(dcmObj, 'UpdateFcn', @myfunction);  % Set the data cursor mode object update
							%   function so it uses updateFcn.m
	% draw Cq position
	Cqx = log(threshold/nInitDnaCopy(CopyRound))/log(E+1);
	hold on, plot([Cqx Cqx], [0 max(slideCounts)*1.2])
	fprintf('\n%g copy, expectation Cq value = %g\nmean Cq value = %g\ndifference between mean value and expectation value: %g‰\n', nInitDnaCopy(CopyRound), Cqx, mean(Cq), abs(mean(Cq) - Cqx) / Cqx * 1000);
	PeakIndex = find(slideCounts(nStart:nEnd) == max(slideCounts(nStart:nEnd)));
	fprintf('%g copy, peak position = %g\n\n', nInitDnaCopy(CopyRound), centers(PeakIndex(1)+nStart-1))
	close(wbar);
	toc
end

function output_txt = myfunction(obj, event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
pos = get(event_obj, 'Position');
output_txt = {['X: ', num2str(pos(1), 6)], ...
              ['Y: ', num2str(pos(2), 6)]};
% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ', num2str(pos(3), 6)];
end
end

function [idx5, idx95] = find90(x, y)
% find region including 90% of data
	area = trapz(x, y);
	if abs(area - 1) > 0.01
		y = y / area;
	end
	cum_prob = cumsum(y);
	cum_prob = cum_prob / max(cum_prob);
	idx5 = find(cum_prob >= 0.05, 1);
	idx95 = find(cum_prob >= 0.95, 1);
end