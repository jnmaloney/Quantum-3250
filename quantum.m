%
%      Bloch states
%      2-Level system
%
%



% Number of atoms to simulate
N = 100;



% Operators
N_a1 = diag(N - (0:N));
N_a2 = diag(0:N);
%N_b = diag(0:N);

%a1_p = diag(sqrt(N + 1 - (1:N)), -1);
%a1_m = diag(sqrt(N + 1 - (1:N)), 1);
%a2_p = diag(sqrt(1:N), -1);
%a2_m = diag(sqrt(1:N), 1);
b_m = diag(sqrt(1:N), 1);
b_p = diag(sqrt(1:N), -1);

f = @(x) (sqrt(x) .* sqrt(N - x + 1));
J_m = diag(f(1:N), 1);
J_p = diag(f(1:N), -1);

J_x = 0.5 * (J_p + J_m);
J_y = 0.5i * (J_p - J_m);
J_z = 0.5 * (N_a1 - N_a2);



% Hamiltonain definition
%H_0 = a1_m * a2_p * b_p;
%H_AB = H_0 + ctranspose(H_0);
%H_m = J_m * b_p + J_p * b_m;
%H_m = J_x;

H_m = b_m + b_p;


% Unitary operator
U = @(t, H) (expm(-1i * H * t));



% Initial state
c0 = zeros(N + 1, 1);
c0(1) = 1;
psi0 = c0;
t = 0.5;
psi_t = U(t, H_m) * psi0;



% Spin-coherent state function (for Qfunc)
alpha = @(theta, phi)  (...
    expm(-1i * phi * J_z) * ...
    expm(-1i * theta * J_y) * ...
    psi0);



% Density matrix
rho_A = psi_t * ctranspose(psi_t);



% Q function
Q = @(theta, phi) ( ...
        ctranspose(alpha(theta, phi)) * ...
        rho_A * ...
        alpha(theta, phi) );

    

% Calculate values
theta_steps = 15;
phi_steps = 30;
q = zeros(theta_steps, phi_steps);
for i = 1:theta_steps
    for j = 1:phi_steps
        
        theta = pi * (i - 1) / (theta_steps - 1);
        phi = 2 * pi * (j - 1) / (phi_steps - 1);

        q(i, j) = abs(Q(theta, phi));

    end
end


% Calculate F_AB

%F_0 = 0.5 * N * (N + 2) - 2 * var(J_z^2);
%F_1 = -var(1i * (J_p - J_m))^2;
%F_2 = -var(J_p^2 + J_m^2);

var = @(A, psi) (ctranspose(psi) * A * psi);

gt_steps = 100;
F_AB = zeros(gt_steps);
gt_space = linspace(0, 0.15, gt_steps);
for i = 1:gt_steps
    gt = gt_space(i);
    psi = U(H_m, gt) * alpha(0, 0);
    F_AB(i) = 4 * (var(J_y^2, psi) - var(J_y, psi)^2);
end

% Plot
figure(1);
semilogy([0, 0.15], [N, N], 'k--');
hold on;
semilogy([0, 0.15], [0.5*N^2, 0.5*N^2], 'k--');
semilogy(gt_space, F_AB, 'r');
axis([0, 0.15, 1e0, 1e4]);

% radius = N / 2;
% theta_matrix = linspace(0, pi, theta_steps);
% phi_matrix = linspace(0, 2*pi, phi_steps);
% [thetamesh, phimesh] = meshgrid(theta_matrix, phi_matrix);
% x = radius * sin(thetamesh) .* cos(phimesh);
% y = radius * sin(thetamesh) .* sin(phimesh);
% z = radius * cos(thetamesh);
% 
% 
% figure(17)
% surf(x, y, z, q')
% shading 'interp'
% view(3)
% %view([1, 0, 0])
% %caxis([0, 2])
% axis on
% axis tight
% axis equal
% %axis([-51, 51, -51, 51, -51, 51])
% grid off;

