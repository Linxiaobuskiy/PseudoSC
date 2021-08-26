
clear all
scRef = readtable('ReferenceSC.csv','ReadRowNames',true);
bulkExp = readtable('Bulk.csv','ReadRowNames',true);
CellLabel = readtable('InfoSC.csv','ReadRowNames',true);
gene_len = readtable('gene_len.txt','ReadRowNames',true);

%% scRNA to Bulk
scBulk = scRNA2Bulk(scRef);
tmpBulk = CellCompositionStatistic1(CellLabel);
scBulk.RealComposition = tmpBulk.composition;
scBulk.RealCellNumbers = tmpBulk.cellnumbers;
clear tmpBulk

%% input data prepare
data = DataPretreatment(scRef,scBulk.matrix);
data.RealComposition = scBulk.RealComposition;
data.RealCellNumbers = scBulk.RealCellNumbers;


%% PseudoSC deconvolution
para.vargenesize = 150;
para.m = 10000;
para.cores = 12;
para.times = 120;
para.mode = 'mosek';

result = PseudoSC(data,para,gene_len);
% cell composition
result = CellCompositionStatistic2(result,CellLabel);

%% NMF deconvolution
para.vargenesize = 150;
para.m = 10000;
para.cores = 12;
para.times = 10;
para.iteration = 500;

result_NMF = nmf_main(data,para,gene_len);
% cell composition
result_NMF = CellCompositionStatistic2(result,CellLabel);