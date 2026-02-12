% 给SimSingleDnaPCRSectRatio.m计算误报率和拒识率的子程序
% 碧云天书 2025年8月28日 7:56:45

function [falseRejectRate, falseAlarmRate] = CalcFrrAndFar(sect, nInitDnaCopy, centersRound, filtedCountsRound)
falseRejectRate = zeros(length(nInitDnaCopy), 1);
falseAlarmCounts = zeros(length(nInitDnaCopy), 1);
falseAlarmRate = zeros(length(nInitDnaCopy), 1);
pos1 = zeros(length(nInitDnaCopy), 1);
pos2 = zeros(length(nInitDnaCopy), 1);
for dnaCopy = 1:length(nInitDnaCopy)
	pos = 0;
	% 查找与多一个模板的曲线相交位置的数据索引号
	while pos < length(centersRound(dnaCopy,:))
		pos = pos + 1;
		if centersRound(dnaCopy, pos) >= sect(dnaCopy, 1)
			break;
		end
	end
	pos1(dnaCopy) = pos;
	% 查找与后一个模板的曲线相交位置的数据索引号
	while pos < length(centersRound(dnaCopy,:))
		pos = pos + 1;
		if centersRound(dnaCopy, pos) >= sect(dnaCopy, 2)
			break;
		end
	end
	pos2(dnaCopy) = pos;

	% 第一个交点前的总数据量
	sum1 = sum(filtedCountsRound(dnaCopy, 1:pos1(dnaCopy)));
	% 第一个交点与第二个交点之间的总数据量
	sum2 = sum(filtedCountsRound(dnaCopy, pos1(dnaCopy)+1:pos2(dnaCopy)));
	% 第二个交点之后的总数据量
	if pos2(dnaCopy) < length(centersRound(dnaCopy,:))
		sum3 = sum(filtedCountsRound(dnaCopy, pos2(dnaCopy)+1:end));
	else
		sum3 = 0;
	end

	% 计算拒识率
	falseRejectRate(dnaCopy) = (sum1 + sum3) / (sum1 + sum2 + sum3);

	% 计算误报率
	falseAlarmCounts(dnaCopy) = 0;
	for ii = 1:length(nInitDnaCopy)					% 循环累积所有dnaCopy以外的模板数扩增后落在dnaCopy有效区间[pos1 pos2]内的数据量
		if ii == dnaCopy			% 自己的曲线不计算
			continue;
		end
		for kk = 1:length(centersRound(ii,:))
			if centersRound(ii, kk) >= sect(dnaCopy, 1) && centersRound(ii, kk) < sect(dnaCopy, 2)
				falseAlarmCounts(dnaCopy) = falseAlarmCounts(dnaCopy) + filtedCountsRound(ii, kk);
			end
		end
	end
	falseAlarmRate(dnaCopy) = falseAlarmCounts(dnaCopy) / (sum1 + sum2 + sum3);
end

% 输出拒识率和误报率
for dnaCopy = 1:length(nInitDnaCopy)
	fprintf('%d copy 误报率 False alarm rate = %f%%\n', dnaCopy, falseAlarmRate(dnaCopy)*100);
end
for dnaCopy = 1:length(nInitDnaCopy)
	fprintf('%d copy 拒识率 False rejection rate = %f%%\n', dnaCopy, falseRejectRate(dnaCopy)*100);
end

end
