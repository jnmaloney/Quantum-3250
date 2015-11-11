function [F_A, F_AB] = qfi(N, M, H, J_x, J_y, J_z, psi0)
% QFI Quantum Fisher Information
%   

    quantum_base;

    % Plot the QFI

    %figure
    %t_steps = 100;
    %t_range = linspace(0, t, t_steps);
    %F_A = zeros(t_steps, 1);
    %F_AB = zeros(t_steps, 1);

    %for i = 1:t_steps

        %t = t_range(i);
        %psi_t = U(t, H) * psi0;

        %rho_AB = psi_t * psi_t';
        %rho_A = eye(N+1).*rho_AB;

        F_AB = 4 * var_d(J_y, rho_AB);
        F_A = 4 * var_d(J_y, rho_A);

    %end

    %semilogy(t_range, F_AB, 'r--');
    %hold on;
    %semilogy(t_range, F_A, 'b');


end

