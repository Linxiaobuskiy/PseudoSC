clear all
%% input data
%{
% scRef: scRNA-seq reference data
scRef = readtable('ReferenceSC.csv','ReadRowNames',true);
% bulkExp: bulk RNA-seq data
bulkExp = readtable('Bulk.csv','ReadRowNames',true);
%}

% Considering the data size limitation in GitHub
% we prepare 'T.mat' to replace scRef and bulkExp
load('T.mat');

% CellLabel: Cell annotations of single cells in scRef
CellLabel = readtable('InfoSC.csv','ReadRowNames',true);
% gene_len: The number of bases of each gene
gene_len = readtable('gene_len.txt','ReadRowNames',true);

%% scRNA to Bulk
scBulk = scRNA2Bulk(scRef); % construct pseudo bulk RNA-seq data from scRNA-seq data
tmpBulk = CellCompositionStatistic1(CellLabel);
scBulk.RealComposition = tmpBulk.composition;
scBulk.RealCellNumbers = tmpBulk.cellnumbers;
clear tmpBulk

%% input data prepare
data = DataPretreatment(scRef,scBulk.matrix);  % construct input data for deconvolution, containing reference matrix and Bulk expression matrix
data.RealComposition = scBulk.RealComposition; % not necessary
data.RealCellNumbers = scBulk.RealCellNumbers; % not necessary


%% PseudoSC deconvolution
para.vargenesize = 150 ;    % parametric genes numbers
para.m = 10000;             % single cell sampling size in each process
para.cores = 12;            % cores used in parallel computing
para.times = 120;           % repeat times
para.mode = 'mosek';        % solver selection

result = PseudoSC(data,para,gene_len); % result.X contains the deconvolution result. X[i,j] means the estimated cell number of cell j in sample i.
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
