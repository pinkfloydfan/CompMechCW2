clear
clc

% constants
EI = 1e2;
EA = 1e5;
q = 1e3;
L = 1;

n = 100;


K = generateGlobalK(n, EI);
F = generateGlobalF(n, q);
x = linspace(0, 1, n);
exact_displacements = q/(24*EI)*x.^2.*(L^2 - 2*L.*x + x.^2);
exact_rotations     = q/(24*EI)*(2*L^2*x - 6*L*x.^2 + 4*x.^3);

rho = zeros(length(F), 1);

L_local = 1/n; 

eps = 100;

while eps > 0.00001 
    
    F_g = evalGlobalF_g(rho, EA, L_local)';
    K_g = evalGlobalK_g(rho, EA, L_local)';
    
    G       = K*rho - F - F_g;
    G_deriv = K + K_g;
    
    %enforcing rho_1, rho_2, rho_n-1, rho_n = 0
    rho_reduced = rho(3:end-2);
    G_reduced = G(3:end-2);
    
    G_deriv_reduced = G_deriv(3:end-2,3:end-2);
    
    new_rho = [0; 0; rho_reduced - G_deriv_reduced\G_reduced ;0 ; 0];
    

    % enforce boundary conditions at each cantilevered end
    %new_rho(1:2,1) = [0;0];
    %new_rho(end-1:end,1) = [0;0];
    
    new_F_g = evalGlobalF_g(new_rho, EA, L_local)';
    new_G   = K*new_rho - F - new_F_g;

   eps = max(abs(new_F_g-F_g));
   
   rho = new_rho;
    
end

fe_displacements = rho(1:2:end-2);
fe_rotations = rho(2:2:end-1); 
%figure 
%loglog(N, errors, '-o')
%grid on
[fe_displacements_linear, fe_rotations_linear] = evalLinearDisplacements(n, EI, q);

figure
hold on
plot(x, fe_displacements_linear, 'DisplayName', 'Linear FE approximation')
plot(x, fe_displacements, 'DisplayName', 'Nonlinear FE approximation')
plot(x, exact_displacements, 'DisplayName', 'Exact solution')
title('Displacements')
grid on 
legend

figure
hold on
plot(x, fe_rotations_linear, 'DisplayName', 'Linear FE approximation')
plot(x, fe_rotations, 'DisplayName', 'Nonlinear FE approximation')
plot(x, exact_rotations, 'DisplayName', 'Exact solution')
title('Rotations')
grid on 
legend

function [fe_displacements_linear, fe_rotations_linear] = evalLinearDisplacements(n, EI, q)

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
    fe_displacements_linear = [0 rho(1:2:end-1)' 0]; 
    fe_rotations_linear     = [0 rho(2:2:end)' 0]; 
end

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

function F_g = evalGlobalF_g(rho, EA, L_local)
    L = length(rho);
    F_g = zeros(1,L);
    for i = 1:2:L-3
        rho_e = rho(i:i+3);
        F_g_e = Fgeom(rho_e, EA, L_local);
        F_g(i:i+3) = F_g(i:i+3) + F_g_e;
    end
end

function K_g = evalGlobalK_g(rho, EA, L_local)
    L = length(rho);
    K_g = zeros(L);
    for i = 1:2:L-3
        rho_e = rho(i:i+3);
        K_g_e = Kgeom(rho_e, EA, L_local);
        K_g(i:i+3,i:i+3) = K_g(i:i+3,i:i+3) + K_g_e;
    end
end
