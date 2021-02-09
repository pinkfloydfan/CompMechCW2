clear
clc

% constants
n = 300; 
EI = 1e2;
EA = 1e5;
q = 1e3;
L = 1;


K = generateGlobalK(n, EI);
F = generateGlobalF(n, q);

% enforce boundary conditions at each cantilevered end
len = length(F);

K_reduced = K(3:len-3,3:len-3);
F_reduced = F(3:len-3, 1);

% calculate displacement vector through simple Ax = b
rho = K_reduced\F_reduced;

% displacements are the odd index terms and rotations are the even index terms
% rho_1 = rho_n = 0 as a part of the boundary conditions
fe_displacements = [0 rho(1:2:end-1)' 0]; 
rotations = [0 rho(2:2:end)' 0]; 

% generate domain
x = linspace(0, 1, n);
exact_displacements = q/(24*EI)*x.^2.*(L^2 - 2*L.*x + x.^2);

figure
hold on
plot(x, fe_displacements, 'DisplayName', 'FE approximation')
plot(x, exact_displacements, 'DisplayName', 'Exact solution')
grid on 
legend


% assemble global stiffness matrix
function K = generateGlobalK(n, EI)
    K = zeros(2*n + 2);
    
    %unit length
    L_local = 1/n; 
    
    K_local = Kmat(EI, L_local);
    
    for i = 1:2:(n*2)
        K(i:i+3,i:i+3) = K(i:i+3,i:i+3) + K_local;
    end
end

% assemble global load vector
function F = generateGlobalF(n, q)
    F = zeros(2*n + 2, 1);
    
    % unit length
    L_local = 1/n; 
    
    for i = 1:2:(n*2)
        F(i:i+3,1) = F(i:i+3,1) + [q*L_local/2 0 q*L_local/2 0]';
    end
end

