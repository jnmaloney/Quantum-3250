function example_discord()
% discord Quantum Discord
%
%
%
%

% ref: PhysRevLett.88.017901.pdf
% Quantum Discord: A Measure of the Quantumness of Correlations
    
    N = 40;
    M = 44;
    
    z_range = linspace(0.0, 0.99, N);
    theta_range = linspace(-pi, pi, M);
    [theta_mesh, z_mesh] = meshgrid(theta_range, z_range);
    q = zeros(N, M);
    
    
    % Parameterised measurement basis 
    Pr0 = @(t) [cos(t); exp(1i)*sin(t)];
    Pr1 = @(t) [exp(-1i)*sin(t); -cos(t)];
    
    % State A
    rho_A = 0.5 * [1 0; 0 1];
    S_A = vn_entropy(rho_A);

    % Calculate all discord values
    for i = 1:N
        for j = 1:M
    
            z = z_range(i);
            theta = theta_range(j);
            
            % Basis
            b0 = Pr0(theta);
            b1 = Pr1(theta);
            
            % Projectors
            P0 = kron(eye(2), b0 * b0');
            P1 = kron(eye(2), b1 * b1');
            
            % Werner state (AB)
            %rho_AB = 0.25 * (1 - z) * eye(4) + ...
            %         0.5 * z * [1 0 0 1]' * [1 0 0 1];
            
            rho_AB = 0.5 * [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1] + ...
                     0.5 * z * [0 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 0];
            
            % Calculate entropies
            m0 = P0 * rho_AB * P0; % / tr_{A,B} (rho_AB * P0)
            m1 = P1 * rho_AB * P1;
            
            Phi = m0 + m1;
            
            S_AB = vn_entropy(rho_AB);
            S_Phi = vn_entropy(Phi);
            
            % Total quantum discord for (z, theta)
            q(i, j) = S_A - S_AB + S_Phi;
                       
        end
    end
    
    figure;
    colormap(spring);
    surf(theta_mesh, z_mesh, q)


end



function S = vn_entropy(p)

    tr = trace(p);
    l = eig(p);
    
    ll = l .* log(l);
    ll(isnan(ll)) = 0;
           
    S = -sum(ll) / tr + log(tr);

    S = real(S);
    
end

