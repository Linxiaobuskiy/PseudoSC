function scBulk = scRNA2Bulk(scRef)

%% matrix
cellname = scRef.Properties.VariableNames;
for i = 1:length(cellname)
    pos = strfind(cellname{i},'_');
    label{1,i} = cellname{i}(1:pos-1);
end

clear pos cellname
list = tabulate(label);
[ncol,~] = size(list);
[nrow,~] = size(scRef);
sz = [nrow,ncol];
for i = 1:ncol
    varTypes{1,i} = 'double';
end
scBulk.matrix = table('Size',sz,'VariableTypes',varTypes,'VariableNames', list(:,1)','RowNames',scRef.Properties.RowNames);
clear sz varTypes varNames nrow

tmp = label;
clear label
[~,label] = ismember(tmp,list(:,1));
clear tmp

matrix = scRef{:,:};

for i = 1:ncol
    flag = ismember(label,i);
    tmp_matrix = matrix(:,flag);
    scBulk.matrix{:,i} = sum(tmp_matrix,2);
end
clear tmp_matrix matrix