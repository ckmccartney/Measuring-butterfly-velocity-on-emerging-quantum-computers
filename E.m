function W = E(G,l,n)
% Creates a brickwall circuit of l layers on n qubits,
    W=1; 
    for j=1:l
        if mod(j,2)== 1
            W=W*LayerA(G(:,:,j),n);
        else
            W=W*LayerB(G(:,:,j),n);
        end
    end
end
