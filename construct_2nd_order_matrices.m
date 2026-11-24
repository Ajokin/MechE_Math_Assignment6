%build the mass and stiffness matrices that describe the 2nd order system.
%M*d2Udt2 + C*dUdt + K*U = B1*uf(t) + B2*duf/dt (eqn 22)
%INPUTS
%string_params: a struct containing the system parameters describing the string
%   string_params.n: number of masses
%   string_params.M: total mass attached to the string
%   string_params.Uf_func: function describing motion of end point
%   string_params.dUfdt_func: time derivative of Uf
%   string_params.Tf: %tension in string
%   string_params.L: %length of string
%   string_params.c: %damping coefficient
%   string_params.dx: %horizontal spacing between masses
%OUTPUTS
%M_mat: the n x n mass (inertia) matrix
%K_mat: the n x n stiffness matrix
%C_mat: (optional) the n x n damping matrix
%B1: (optional) n x 1 column, how uf(t) pushes on the masses
%B2: (optional) n x 1 column, how duf/dt pushes on the masses
function [M_mat,K_mat,C_mat,B1,B2] = construct_2nd_order_matrices(string_params)
    n = string_params.n;
    M = string_params.M;
    Tf = string_params.Tf;
    c = string_params.c;
    dx = string_params.dx;

    %construct the nxn discrete laplacian matrix
    I_n = eye(n); % build the nxn identity matrix

    %-2 down the diagonal, then shift the identity right and left to get
    %the 1's right above and right below the diagonal (eqn 16)
    my_Laplacian = -2*I_n + circshift(I_n,[0,1]) + circshift(I_n,[0,-1]);

    %doing this as subtracting deals with the weird edge cases
    %that occur when n=1 and n=2
    my_Laplacian(1,end) = my_Laplacian(1,end)-1; %delete unwanted 1 in top right corner
    my_Laplacian(end,1) = my_Laplacian(end,1)-1; %delete unwanted 1 in bottom left corner

    %nth cartesian basis vector, only the last mass is touching the actuator
    e_n = I_n(:,n);

    %eqn 23
    M_mat = (M/n)*I_n;
    K_mat = -(Tf/dx)*my_Laplacian;
    C_mat = -(c/dx)*my_Laplacian;
    B1 = (Tf/dx)*e_n;
    B2 = (c/dx)*e_n;
end
