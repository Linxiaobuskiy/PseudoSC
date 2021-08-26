function data = nmf_main(data,para,gene_len)

% m, sampling size
% K, iteration times

data.X = zeros(data.n_bulk,data.n_cells);
A = 1:data.n_cells;

%% normalization
E = data.matrix.bulk(1:para.vargenesize,:);
G = data.matrix.reference(1:para.vargenesize,:);

for i = 1:para.vargenesize
    E(i,:) = E(i,:)/gene_len{data.ParametricGenes{i},'Length'};
    G(i,:) = G(i,:)/gene_len{data.ParametricGenes{i},'Length'};
end

%% NMF
V = E';
clear E;
Ref = G';
clear G;

parfor i = 1:para.times
    random_num(i,:) = A(randperm(numel(A),para.m));
    random_num(i,:) = sort(random_num(i,:));
    H = Ref(random_num(i,:),:);
    W{i,1} = nmf_sc(V,H,para.iteration);
    fprintf('Processing for repeat times %i ;\n',i);
end

for i = 1:para.times
    data.X(:,random_num(i,:)) = data.X(:,random_num(i,:)) + W{i,1}/para.times;
end
clear W;

for i = 1:data.n_bulk
    data.distance(i) = (data.X(i,:)*Ref-V(i,:))*(data.X(i,:)*Ref-V(i,:))';
end
