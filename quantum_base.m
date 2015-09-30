% Quantum Mechanics function set


% Unitary operator
U = @(t, H) (expm(-1i * H * t));


% Variance function
%var = @(A, psi) ((psi' * (A^2) * psi) - (psi' * A * psi)^2);


% Variance by density operator
var_d = @(A, rho) (trace(rho * A^2) - trace(rho * A)^2);


% Spin-coherent state function (for Qfunc)
alpha = @(theta, phi)  (...
    expm(-1i * phi * J_z) * ...
    expm(-1i * theta * J_y) * ...
    psi0);


% Q function
Q = @(theta, phi) ( ...
        alpha(theta, phi)' * ...
        rho_A * ...
        alpha(theta, phi) );

    
% von Neumann entropy
vn_entropy = @(rho) -trace(rho * logm(rho));


% basis projector
%projector = @(i) sparse(i, i, 1, N + 1, N + 1); 
el = @(b, i) b(i);
b0 = @(theta) [cos(theta), sin(theta) * exp(1i)];
b1 = @(theta) [sin(theta) * exp(1i), cos(theta)];
b = @(i, theta) [el(b0(theta), i); el(b1(theta), i)];
projector = @(i, theta) b(i, theta) * b(i, theta)';


% Probability of outcome i
P = @(i, theta, rho) trace(projector(i, theta) * rho);


% Projected state of rho for P = P_i
projected_state = @(i, theta, rho) projector(i, theta) * rho * projector(i, theta) / trace(projector(i, theta) * rho);


% Projector expression ?
f = @(i, theta, rho) P(i, theta, rho) * vn_entropy(projected_state(i, theta, rho));
projected_rho = @(theta, rho) sum(f(1:N, theta, rho));
% P(i) *
% vn_entropy(projected_state(projector(i), rho_A))
%proj = @(A, B, AB) ();



% Discord I = H(S:A)
%discord_I = @(S, A) (vn_entropy(S) + vn_entropy(A) - vn_entropy(S, A));
%discord_I = @(A, B, AB) (vn_entropy(A) + vn_entropy(B) - vn_entropy(AB));


% Discord J = H(S|A)
%discord_J = @(S, A) (vn_entropy(S) + vn_entropy(A) - vn_entropy(S, A));


% Discord
%discord = @(S, A) (discord_I(S, A) - discord_J(S, A));
discord = @(A, AB, theta) (vn_entropy(A) - vn_entropy(AB) + vn_entropy(projected_rho(theta, AB)));

