%
%      Bloch states
%      2-Level system
%
%



quantum_base;



% Number of atoms to simulate
N = 1e2 + 1;
M = 0;

% Evolution Time
t = 0.40;



% Initial state
psi0 = [1; zeros(N, 1)];
%psi0 = 1/sqrt(N+1) * ones(N + 1, 1);
%psi0 = 1/sqrt(2) * [1; zeros(N, 1); 1; zeros(N, 1)];

% Operators
N_a1 = diag(N - (0:N));
N_a2 = diag(0:N);
%N_b = diag(0:N);

a1_p = diag(sqrt(N + 1 - (1:N)), -1);
a1_m = diag(sqrt(N + 1 - (1:N)), 1);
a2_p = diag(sqrt(1:N), -1);
a2_m = diag(sqrt(1:N), 1);

f = @(x) (sqrt(x + M));
b_m = diag(f(1:N), 1);
b_p = diag(f(1:N), -1);

f = @(x) (sqrt(x) .* sqrt(N - x + 1));
J_p = diag(f(1:N), 1);
J_m = diag(f(1:N), -1);

J_x = 0.5 * (J_p + J_m);
J_y = 0.5i * (J_p - J_m);
J_z = 0.5 * (N_a1 - N_a2);



% Hamiltonain definition
H_m = J_m .* b_p + J_p .* b_m;



% Psi0 at time t
psi_t = U(t, H_m) * psi0;



% Density matrix
rho_AB = psi_t * psi_t';
rho_A = eye(N+1).*rho_AB;
meanJz = trace(J_z*rho_A);
varJy = trace(J_y*J_y*rho_A) - trace(J_y*rho_A).^2;



% Calculate the QFI of the state after evolution
H = H_m;

% figure
% hold on;
% 
t_steps = 100;
t_range = linspace(0, t, t_steps);
F_A = zeros(t_steps, 1);
F_AB = zeros(t_steps, 1);
D = zeros(t_steps, 1);

for i = 1:t_steps

    t = t_range(i);
    psi_t = U(t, H) * psi0;

    rho_AB = psi_t * psi_t';
    rho_A = eye(N+1).*rho_AB;

    F_AB(i) = 4 * var_d(J_y, rho_AB);
    F_A(i) = 4 * var_d(J_y, rho_A);

    % Calculate the discord of the state if possible
    D(i) = qd(N, J_x, J_y, rho_AB, rho_A);

end

%title(sprintf('t = %f', t));

%semilogy(t_range, F_AB / N, 'r--');
%semilogy(t_range, F_A / N, 'k');
%semilogy(t_range, D, 'b');

% QD v QFI
figure
hold on
title(sprintf('t = %f', t));

% Create inset discord / bloch
%axes('Parent', hf1, 'Position',[0.20 0.65 0.25 0.25]);
%box on
subplot(1, 3, 3);
colormap(spring)

bloch_inset(N, M, t, H, J_x, J_y, J_z, psi0, rho_A);

subplot(1, 3, 2);
hold on;

semilogy(t_range, F_AB / N, 'r--');
semilogy(t_range, F_A / N, 'k');
semilogy(t_range, D, 'b');

subplot(1, 3, 1);
    axis on
    axis tight
    %axis equal
    grid off;
    xlabel('F_AB')
    zlabel('D')
    set(gca, 'XTickLabel', '');
    set(gca, 'YTickLabel', '');
    set(gca, 'ZTickLabel', '');
plot(F_AB / N, D);


