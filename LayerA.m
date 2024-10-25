function L1 = LayerA(g, n)
    %% Layer1(g,n) = g*g*g..*g (n even) g*g*g..g*I (n odd)
    % where * denotes kroneker tensor product
    L1=1;
    for j = 2:2:n
        L1 = kron(L1, g);
    end
    if mod(n,2)==1
        L1 = kron(L1,eye(2));
    end
end