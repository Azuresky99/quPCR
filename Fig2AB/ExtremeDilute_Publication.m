%   Smooth for ExtremeDiluteR.xlsx
%		2025-7-22

% Parameters
histNbin = 1800;								% bin value of histgram
profileWidth = 1;								% smooth width

%% import extreme dilute qPCR data
filename = "ExtremeDiluteR.xlsx";
opts = spreadsheetImportOptions("NumVariables", 6);
Sheet1Name = "ExtreDial";
opts.Sheet = Sheet1Name;
opts.DataRange = "A2:F700";
opts.VariableNames = ["Well", "Despription", "FluoChn", "TargetGene", "PrimeNum", "Cq"];
opts.SelectedVariableNames = "Cq";
opts.VariableTypes = ["string", "string", "string", "string", "string", "double"];
opts.MissingRule = "omitrow";
opts = setvaropts(opts, ["Well", "Despription", "FluoChn", "TargetGene", "PrimeNum"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Well", "Despription", "FluoChn", "TargetGene", "PrimeNum"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "Cq", "TreatAsMissing", '');
SingleChainCq = readtable(filename, opts, "UseExcel", false);
SingleChainCq = table2array(SingleChainCq);
clear opts

%% Process data
% show Cq distribution
figure, hist(SingleChainCq, histNbin)
[counts1,centers1] = hist(SingleChainCq, histNbin); title(sprintf('Cq distribution, %d replication wells in total', length(SingleChainCq)));
% smooth by slide window. Window width set at 0.15 cycles
slideWin = round(0.15 / (centers1(2) - centers1(1)));
slideWin = floor(slideWin/2) * 2 + 1;			% make sure window is odd
% add zero data before and after counts1
extCount = slideWin;
extZero = zeros(1, extCount);
countsExt = [extZero, counts1, extZero];
nStart = floor(extCount / 2) + 1;
nEnd = length(countsExt) - floor(extCount / 2);
% add width of half window before and after centers1
dCenter = centers1(2)-centers1(1);
centersExtHead = (centers1(1)-ceil(extCount/2)*dCenter) : dCenter : (centers1(1)-dCenter);
centersExtTail = (centers1(end)+dCenter) : dCenter : (centers1(end)+ceil(extCount/2)*dCenter);
centersExt = [centersExtHead, centers1, centersExtTail];
slideCounts1 = zeros(size(counts1));
halfWidth = floor(slideWin / 2);
% create Gaussian distribution curve
ntemplate = normpdf(1:slideWin, nStart, slideWin/5);
% Gaussian smooth
for ii = nStart:nEnd
	slideCounts1(ii-floor(extCount/2)) = sum(countsExt((ii-halfWidth):(ii+halfWidth)) .* ntemplate);
end
% Plot figures
figure, bar(centers1, counts1, 1), hold on, plot(centersExt, slideCounts1, 'r', 'linewidth', profileWidth), title('Cq distribution + Cq distribution after smooth'), legend("Cq", "smooth Cq"), xlim([35,50])
figure, bar(centersExt, slideCounts1, 1), title('Cq distribution after smooth'), xlim([35,50])
% set precision cursor
dcmObj = datacursormode;  % Turn on data cursors and return the data cursor mode object
set(dcmObj, 'UpdateFcn', @myfunction);  % Set the data cursor mode object update function so it uses updateFcn.m

%% precision cursor
function output_txt = myfunction(obj, event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
pos = get(event_obj, 'Position');
output_txt = {['X: ', num2str(pos(1), 6)], ['Y: ', num2str(pos(2), 6)]};
% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ', num2str(pos(3), 6)];
end
end