function D = qd(N, J_x, J_y, rho_AB, rho_A)
% Quantum discord
%   
    quantum_base;

    theta = 0;
    phi = 0;
    
    %theta_steps = 20;
    %phi_steps = 20;
        
    %phi_range = linspace(-pi, pi, theta_steps);
    %theta_range = linspace(-pi, pi, theta_steps);
    %[theta_mesh, phi_mesh] = meshgrid(theta_range, phi_range);
    %q = zeros(phi_steps, theta_steps);
    
    
    % Parameterised measurement basis 
    Pr0 = @(k, theta, phi) circshift(eye(N + 1, 1), k, 1);
    
    % State A
    rho_A = 0.5 * [1 0; 0 1];
    S_A = vn_entropy(rho_A);

    % Calculate all discord values
    %for i = 1:phi_steps
    %    for j = 1:theta_steps
    
    %        phi = phi_range(i);
    %        theta = theta_range(j);
            
            Phi = 0;
            for k = 1:N
                b0 = Pr0(k, theta, phi);
                P0 = b0 * b0';
                %size(rho_AB)
                %size(P0)
                m0 = P0 * rho_AB * P0;
                Phi = Phi + m0;
            end
            
            S_AB = vn_entropy(rho_AB);
            S_Phi = vn_entropy(Phi);
            
            % Total quantum discord for (z, theta)
            %q(i, j) = S_A - S_AB + S_Phi;
         
    %    end
    %end
    
    %figure;
    %colormap(spring);
    %surf(theta_mesh, phi_mesh, q)
    
    D = S_A - S_AB + S_Phi;
    
end



function V = vn_entropy(p)

    tr = trace(p);
    l = eig(p);
    
    ll = l .* log(l);
    ll(isnan(ll)) = 0;
           
    V = -sum(ll) / tr + log(tr);

    V = real(V);
    
end

