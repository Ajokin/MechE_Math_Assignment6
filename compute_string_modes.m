%Builds the mass and stiffness matrices for the string and gets the mode
%shapes and resonant frequencies from the generalized eigenvalue problem
%K*Ur = lambda*M*Ur, where lambda = omega_r^2 (eqn 26)
%(no damping and no forcing, so both ends of the string are held at 0)
%INPUTS:
%string_params: a struct containing the system parameters describing the string
%OUTPUTS:
%Ur_mat: n x n matrix where column i is the ith mode shape
%omega_n: column of resonant frequencies, sorted lowest to highest
function [Ur_mat,omega_n] = compute_string_modes(string_params)
    [M_mat,K_mat] = construct_2nd_order_matrices(string_params);

    %Use MATLAB to solve the generalized eigenvalue problem
    [Ur_mat,lambda_mat] = eig(K_mat,M_mat);

    %eigenvalues are omega^2, so take the square root to get the
    %resonant frequencies
    omega_n = sqrt(diag(lambda_mat));

    %put them in order lowest to highest, just in case eig doesn't
    [omega_n,sort_idx] = sort(omega_n);
    Ur_mat = Ur_mat(:,sort_idx);

    %a mode shape is only defined up to a scale factor, so scale each one
    %so its biggest entry is exactly +1 (makes the plots easier to compare)
    for i = 1:string_params.n
        [~,max_idx] = max(abs(Ur_mat(:,i)));
        Ur_mat(:,i) = Ur_mat(:,i)/Ur_mat(max_idx,i);
    end
end
