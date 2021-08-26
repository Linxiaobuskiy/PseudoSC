function W = nmf_sc(V,H,K)
    [m,~]=size(V);  
    [R,~]=size(H);
    W = abs(rand(m,R));
    a = V*H';
    b = H*H';
    for i = 1:K
        %fprintf('Iteration %i; \n',i);
        %W = W .* (V*H') ./ (W*(H*H'));
        W = W .* a ./ (W*b);
        %loss(i) = sum(sum((V-W*H).*(V-W*H)));
    end
end  