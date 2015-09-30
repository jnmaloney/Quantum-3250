%
%      Bloch states
%      2-Level system
%
%



quantum_base;



% Number of atoms to simulate
N = 1;
M = 0;

% Evolution Time
t = 0.15;



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


% ...
qfi(N, M, t, H_m, J_x, J_y, J_z, psi0, rho_A);


