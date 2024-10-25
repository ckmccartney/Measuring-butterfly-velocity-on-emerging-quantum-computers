function L1 = LayerB(g, n)
    %% Layer1(g,n) = I*g*g*g..*g*I (n even) I*g*g*g..g (n odd)
    % where * denotes kroneker tensor product

    L1 = eye(2);
    for j = 3:2:n
        L1 = kron(L1, g);
    end
    if mod(n,2)==0
        L1 = kron(L1,eye(2));
    end
end