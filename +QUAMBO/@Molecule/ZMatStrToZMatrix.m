function zmatrix = ZMatStrToZMatrix(zMatStr)
zMatStr = regexprep(zMatStr, '[,;\n\r\n]', ';');
indDelim = regexp(zMatStr, ';');

if(length(indDelim)<2)
    zMatStr = [zMatStr(1:indDelim(1)-1), ' 0 0 0 0 0 0;', zMatStr(indDelim(1)+1:end)];
elseif(length(indDelim)<3)
    zMatStr = [zMatStr(1:indDelim(1)-1), ' 0 0 0 0 0 0;', ...
        zMatStr(indDelim(1)+1:indDelim(2)-1), ' 0 0 0 0;', zMatStr(indDelim(2)+1:end)];
else
    zMatStr = [zMatStr(1:indDelim(1)-1), ' 0 0 0 0 0 0;', ...
        zMatStr(indDelim(1)+1:indDelim(2)-1), ' 0 0 0 0;', ...
        zMatStr(indDelim(2)+1:indDelim(3)-1), ' 0 0;', zMatStr(indDelim(3)+1:end)];
end

periodicTable = Molecule.InitializePeriodicTable();
for number = 1:length(periodicTable)
    symbol = periodicTable{number};
    zMatStr = regexprep(zMatStr, [symbol,' '], [num2str(number),' ']);
end
zmatrix = str2num(zMatStr);
end