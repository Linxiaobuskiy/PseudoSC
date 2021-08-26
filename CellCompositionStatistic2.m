function data = CellComposistion(data,CellLabel)

celltype = tabulate(CellLabel{:,:});
celltype = sortrows(celltype,1);

[nrow,~] = size(celltype);
ncol = data.n_bulk;
sz = [nrow,ncol];
for i = 1:ncol
    varTypes{1,i} = 'double';
end
data.EstiCellNumbers = table('Size',sz,'VariableTypes',varTypes,'VariableNames',data.BulkSampleNames,'RowNames',celltype(:,1)');
data.EstiComposition = table('Size',sz,'VariableTypes',varTypes,'VariableNames',data.BulkSampleNames,'RowNames',celltype(:,1)');

for i = 1:nrow
    index{i} = matches(CellLabel{:,:},celltype{i,1});
end



for i = 1:data.n_bulk
    tmp = data.X(i,:);
    counter = 0;
    for j = 1:nrow
        data.EstiCellNumbers{j,i} = sum(tmp(index{j}));
        counter = counter+data.EstiCellNumbers{j,i};
    end
    for j = 1:nrow
        data.EstiComposition{j,i} = data.EstiCellNumbers{j,i}/counter;
    end
end