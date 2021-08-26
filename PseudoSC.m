function output = PseudoSC(data,para,gene_len)

% para: parameters
% accuracy: High or Normal
m = para.m;
times = para.times;
cores = para.cores;
mode = para.mode;

%if cores>1
%    parpool(cores)
%end


A = 1:data.n_cells;

data.X = zeros(data.n_bulk,data.n_cells);
data.estiCellNumber = zeros(1,data.n_bulk);
data.RepeatTimes = zeros(1,data.n_bulk);
data.SolvedTimes = zeros(1,data.n_bulk);

%% %% normalization
Bulk = data.matrix.bulk(1:para.vargenesize,:);
Ref = data.matrix.reference(1:para.vargenesize,:);

for i = 1:para.vargenesize
    Bulk(i,:) = Bulk(i,:)/gene_len{data.ParametricGenes{i},'Length'};
    Ref(i,:) = Ref(i,:)/gene_len{data.ParametricGenes{i},'Length'};
end

%% cvx
for i = 1:data.n_bulk
    
    E = Bulk(:,i)';
    tmp.X = zeros(1,data.n_cells);
    
    if cores>1
        
        while data.SolvedTimes(i) < times
            parfor j = 1:times
                %fprintf('processing sample %d for times %d \n',i,j);
                random_num(j,:) = A(randperm(numel(A),m));
                random_num(j,:) = sort(random_num(j,:));
                G = Ref(:,random_num(j,:));
                %[tmp(j),cvx_status{j},cvx_optval(j),c(j)] = core_cvx(E,G,1);
                [cvx_X{j},cvx_status{j},~,c(j)] = cvxDC(E,G,0,mode);
                fprintf('processing sample %d for times %d ; cvx status: %s ;\n',i,data.RepeatTimes(i)+j,cvx_status{j});
            end
            
            data.RepeatTimes(i) = data.RepeatTimes(i)+times;
            
            for j = 1:times
                if (strcmp(cvx_status{j},'Solved')==1)||(strcmp(cvx_status{j},'Inaccurate/Solved')==1)
                    data.SolvedTimes(i) = data.SolvedTimes(i) + 1;
                    tmp.X(data.SolvedTimes(i),random_num(j,:)) = cvx_X{j};
                    tmp.estiCellNumber(data.SolvedTimes(i)) = c(j);
                    XX = tmp.X(data.SolvedTimes(i),:);
                    tmp.distance(data.SolvedTimes(i)) = (XX*Ref'-E)*(XX*Ref'-E)';
                end
            end
        end
        
    else    % cores == 1
        while data.SolvedTimes(i) < times
            data.RepeatTimes(i) = data.RepeatTimes(i)+1;
            j = data.RepeatTimes(i);
            random_num(j,:) = A(randperm(numel(A),m));
            random_num(j,:) = sort(random_num(j,:));
            G = Ref(:,random_num(j,:));
            %[tmp(j),cvx_status{j},cvx_optval(j),c(j)] = core_cvx(E,G,1);
            [cvx_X{j},cvx_status{j},~,c(j)] = core_cvx(E,G,0,mode);
            fprintf('processing sample %d for times %d ; cvx status: %s ;\n',i,j,cvx_status{j});
            
            if (strcmp(cvx_status{j},'Solved')==1)||(strcmp(cvx_status{j},'Inaccurate/Solved')==1)
                data.SolvedTimes(i) = data.SolvedTimes(i) + 1;
                tmp.X(data.SolvedTimes(i),random_num(j,:)) = cvx_X{j};
                tmp.estiCellNumber(data.SolvedTimes(i)) = c(j);
                XX = tmp.X(data.SolvedTimes(i),:);
                tmp.distance(data.SolvedTimes(i)) = (XX*Ref'-E)*(XX*Ref'-E)';
            end
            
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for j = 1:times
        data.X(i,:) = data.X(i,:) + tmp.X(j,:);
    end
    
    data.X(i,:) = data.X(i,:)/times;
    data.estiCellNumber(i) = sum(data.X(i,:));
    XX = data.X(i,:);
    data.distance(i) = (XX*Ref'-E)*(XX*Ref'-E)';
end

data.para = para;
output = data;