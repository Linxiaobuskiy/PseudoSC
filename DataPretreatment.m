function output = DataPretreatment(scRef,bulkExp)

data.matrix.reference = scRef{:,:};
data.matrix.bulk = bulkExp{:,:};
[~,data.n_cells] = size(data.matrix.reference);
[~,data.n_bulk] = size(data.matrix.bulk);

data.ParametricGenes = bulkExp.Properties.RowNames;
data.BulkSampleNames = bulkExp.Properties.VariableNames;

output = data;