%
%      Bloch states
%      2-Level system
%
%



% Number of atoms to simulate
N = 20;



% Operators
N_a1 = diag(N - (0:N));
N_a2 = diag(0:N);
%N_b = diag(0:N);

a1_p = diag(sqrt(N + 1 - (1:N)), -1);
a1_m = diag(sqrt(N + 1 - (1:N)), 1);
a2_p = diag(sqrt(1:N), -1);
a2_m = diag(sqrt(1:N), 1);
b_p = diag(sqrt(1:N), -1);
b_m = diag(sqrt(1:N), 1);

f = @(x) (sqrt(x) .* sqrt(N - x + 1));
J_m = diag(f(1:N), 1);
J_p = diag(f(1:N), -1);

J_x = 0.5 * (J_p + J_m);
J_y = 0.5i * (J_p - J_m);
J_z = 0.5 * (N_a1 - N_a2);


    
% Initial state
c0 = zeros(N + 1, 1);
c0(1) = 1;

psi0 = c0;



% State function
alpha = @(theta, phi)  (...
    expm(-1i * theta * J_y) * ...
    expm(-1i * phi * J_z) * ...
    psi0);



% Density matrix
rho_A = psi0 * psi0';



% Q function
Q = @(theta, phi) ( ...
        ctranspose(alpha(theta, phi)) * ...
        rho_A * ...
        alpha(theta, phi) );



% Hamiltonain definition
H_0 = a1_m * a2_p * b_p;
H_AB = H_0 + ctranspose(H_0);

H_m = J_m * b_p + J_p * b_m;



% Unitary operator
U = @(t, H) (expm(-1i * H * t));



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



% Plot
radius = N / 2;
theta_matrix = linspace(0, pi, theta_steps);
phi_matrix = linspace(0, 2*pi, phi_steps);
[thetamesh, phimesh] = meshgrid(theta_matrix, phi_matrix);
x = radius * sin(thetamesh) .* cos(phimesh);
y = radius * sin(thetamesh) .* sin(phimesh);
z = radius * cos(thetamesh);


figure(17)
surf(x, y, z, q')
shading 'interp'
view(3)
%view([1, 0, 0])
%caxis([0, 2])
axis on
axis tight
axis equal
%axis([-51, 51, -51, 51, -51, 51])
grid off;

