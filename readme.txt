Introduction:

Tool for bulk RNA-seq data deconvolution.

Environment:
1. MATLAB, version R2020a suggested.
2. CVX: Matlab Software for Disciplined Convex Programming. (http://cvxr.com/cvx/)
3. MOSEK solver. (http://cvxr.com/cvx/doc/mosek.html)
4. Hardware: RAM 128GB, CPU i9-7920X(12 cores).

Input data:
1. Single-cell expression reference matrix.
2. Bulk expression matrix.
3. Cell label annotation file.
4. Parameters (including parametric gene numbers, sampling size, parallel core numbers, sampling times, mode).
5. Exon length annotation file.

functions:
1. scRNA2Bulk: generate pesudo bulk expression from single cell expression.
2. CellCompositionStatistic1: caculate the cell compostition of reference expression.
2. CellCompositionStatistic2: caculate the cell compostition of ouput pseudo single cell expression.
3. DataPretreatment: Pretreat the input data for deconvolution.
4. PseudoSC: main function of deconvolution.
5. cvxDC: The module functions of main function.
6. nmf_main: NMF method for deconvolution.
