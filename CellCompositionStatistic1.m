function scBulk = CellCompositionStatistic(CellLabel)

%% composition
cellname = CellLabel.Properties.RowNames;
for i = 1:length(cellname)
    pos = strfind(cellname{i},'_');
    label{1,i} = cellname{i}(1:pos-1);
end
clear pos cellname
list = tabulate(label);

celltypes = tabulate(CellLabel{:,1});
celltypes = sortrows(celltypes,1);
[nrow,~] = size(celltypes);
[ncol,~] = size(list);

sz = [nrow,ncol];
for i = 1:ncol
    varTypes{1,i} = 'double';
end
scBulk.cellnumbers = table('Size',sz,'VariableTypes',varTypes,'VariableNames', list(:,1)','RowNames',celltypes(:,1)');
scBulk.composition = table('Size',sz,'VariableTypes',varTypes,'VariableNames', list(:,1)','RowNames',celltypes(:,1)');

tmp = label;
clear label
[~,label] = ismember(tmp,list(:,1));
clear tmp

for i = 1:ncol
    flag1 = ismember(label,i);
    tmp_label = CellLabel(flag1,1);
    for j = 1:nrow
        flag2 = ismember(tmp_label{:,:},celltypes(j,1));
        scBulk.cellnumbers{j,i} = sum(flag2);
        scBulk.composition{j,i} = scBulk.cellnumbers{j,i}/list{i,2};
    end
end

