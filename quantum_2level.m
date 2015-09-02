%
%      Expectation values of a
%      single 2-Level system
%
%


% Operator formulation
a_p = [0 0; 1 0];
a_m = [0 1; 0 0];
b_p = [0 0; 1 0];
b_m = [0 1; 0 0];

% Tensor product operator
H_m = kron(a_p, b_m) + kron(a_m, b_p);

% Unitary evolution operator
U = @(t, H) (expm(-1i * H * t));

% Number operators
N_a = kron(eye(2), a_p * a_m);
N_b = kron(b_p * b_m, eye(2));

% State
a0 = [0; 1];
b0 = [1; 0];
x0 = kron(a0, b0);

% Create a series of evolved states
gt_steps = 50;
t_space = linspace(0, 3.0, gt_steps);
Na = zeros(gt_steps, 1);
Nb = zeros(gt_steps, 1);
for i = 1:gt_steps
    t = t_space(i);
    x = U(t, H_m) * x0;
    Na(i) = ctranspose(x) * N_a * x;
    Nb(i) = ctranspose(x) * N_b * x;
end

figure(1);
hold on;
plot(t_space, Na, 'r');
plot(t_space, Nb, 'b');
