clear
clc

digits(100)
% Let,
I = [1 0; 0 1];
X = [0 1; 1 0];
Y = [0 -1i; 1i 0];
Z = [1 0; 0 -1];

% and define the Kiteav Chain on n qubits, so(n) ...
xx = kron(X,X);
XX = @(j,n) kron( eye(2^(j-1)), kron(xx, eye(2^(n-j-1))) );  

yy = kron(Y,Y);
YY = @(j,n) kron( eye(2^(j-1)), kron(yy, eye(2^(n-j-1))) );  

zi = kron(Z,I);
ZI = @(j,n) kron( eye(2^(j-1)), kron(zi, eye(2^(n-j-1))) );  

% ... for n qubits,
n=3;

r = 0;
h = 0;
% .. define H = Σ_j XY(j,n), as 
H = zeros(size(XX(1,n)));
ar = randn(n-1,1);
for j = 1:n-1
    H = H + (1+r)*XX(j,n)/2 + (1-r)*YY(j,n)/2 +  h*ZI(j,n);
end

H = H+ h*kron(eye(2^(n-1)), Z);


%We would like to find a nice circuit decposition to simulate this unitary with ... 

% ... m SU(2) layers
m = 3;

% define the product manifold of m copies of U(2) 
mfld = unitaryfactory(4, m);

% Define the encoding map, E: mfld -> Circuit as E(G,m,n) \\-> E.m

x = zeros(4, 4, m, 11);
for i = 0:10
    t = i/2;
    
    %Let the true unitary be 
    U = expm(-1i*H*t);
    % sequence = randi([1 2], 11, 1);
    
    % Setting up the RTR-Prototype
    problem.M = mfld;
    problem.cost = @(X) - real(trace(U'*E(X, m, n))/2^n);
    
    %Let the initial circuit encoding be x0 
    for j=1:m
        x0(:,:,j) = expm(-1i*(xx+yy)*t/m);
    end
            
        
    % Using automatic differentiation package-method to compute Gradients and Hessians
    problem = manoptAD(problem);
    % checkgradient(problem); 
    
    % Set options for trustregions solver
    options.tolcost = - 0.999; % set tolerance for cost function to 99% accuracy
    options.maxiter = 100; % set maximum number of iterations (optional)
    
    
    % Using manopt package to solver for RTR 
    [x(:, :, :, i+1), xcost, info, options] = trustregions(problem, x0, options);
    disp(t)
end

save('Hamiltonians3Param1.mat', 'x')

% display(norm(E(x, m, n)-U))