%
%      Bloch states
%      2-Level system
%
%



% Number of atoms to simulate
N = 10;
t = 0.04;
M = 2;




% Initial state
c0 = zeros(N + 1, 1);
c0(1) = 1;
psi0 = c0;



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



% Unitary operator
U = @(t, H) (expm(-1i * H * t));



% Spin-coherent state function (for Qfunc)
alpha = @(theta, phi)  (...
    expm(-1i * phi * J_z) * ...
    expm(-1i * theta * J_y) * ...
    psi0);



% Psi0 at time t
psi_t = U(t, H_m) * psi0;



% Density matrix
rho_AB = psi_t * psi_t';
rho_A = eye(N+1).*rho_AB;
meanJz = trace(J_z*rho_A);
varJy = trace(J_y*J_y*rho_A) - trace(J_y*rho_A).^2;

 

% Q function
Q = @(theta, phi) ( ...
        alpha(theta, phi)' * ...
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



% Variance function
var = @(A, psi) (psi' * A * psi);



% Bloch space
theta_matrix = linspace(0, pi, theta_steps);
phi_matrix = linspace(0, 2*pi, phi_steps);
[thetamesh, phimesh] = meshgrid(theta_matrix, phi_matrix);



% Calculate eigenvalues of J_y for P(J_y)
[V, D] = eig(J_y);
P_y = zeros(N + 1, 1);

for i = 1:N + 1
    % ith eigenvector
    v = V(:,i);
    
    P_y(i) = v' * rho_A * v;
end




figure
colormap(spring);
%title(sprintf('Husimi-Q function and J_y projection for N=%l atoms, M=%i photons and gt=%i', N, M, t));



subplot(2, 1, 1)


radius = N / 2;
x = radius * sin(thetamesh) .* cos(phimesh);
y = radius * sin(thetamesh) .* sin(phimesh);
z = radius * cos(thetamesh);

surf(x, y, z, q')
shading 'interp'
view(3)
axis on
axis tight
axis equal
grid off;
     
     

subplot(2, 1, 2)
    
bar(-N/2:N/2, P_y, 'k');
xl = xlabel('J_y');
yl = ylabel('P(J_y)');
set([xl yl], 'interpreter', 'tex');

     
 
 % Plot the QFI
 
 figure
 t_steps = 100;
 t_range = linspace(0, t, t_steps);
 F_A = zeros(t_steps, 1);
 F_AB = zeros(t_steps, 1);
 
 for i = 1:t_steps
 
 t = t_range(i);
 psi_t = U(t, H_m) * psi0;
 
 rho_AB = psi_t * psi_t';
 rho_A = eye(N+1).*rho_AB;
 
 F_AB(i) = 4 * var(J_y, rho_AB);
 F_A(i) = 4 * var(J_y, rho_A);
 
 end
 
 semilogy(t_range, F_AB, 'r--');
 hold on;
 semilogy(t_range, F_A, 'b');
