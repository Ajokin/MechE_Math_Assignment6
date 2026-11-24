%Rate function encoding the equations of motion for the vibrating string
%(n masses on a taut string, left end fixed to the wall and the right end
%moved by the actuator)
%INPUTS
%t: current time
%V: system state. V = [U;dUdt] where
%   U and dUdt are n x 1 column vectors
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
%dVdt = [dUdt;d2Udt2]: the time derivative of the state vector
function dVdt = string_rate_func01(t,V,string_params)
    n = string_params.n; %number of masses
    M = string_params.M; %total mass attached to the string
    Uf_func = string_params.Uf_func; %function describing motion of end point
    dUfdt_func = string_params.dUfdt_func; %time derivative of Uf
    Tf = string_params.Tf; %tension in string
    L = string_params.L; %length of string
    c = string_params.c; %damping coefficient
    dx = string_params.dx; %horizontal spacing between masses

    %unpack state variable
    U = V(1:n);
    dUdt = V((n+1):(2*n));

    Uf = Uf_func(t);
    dUfdt = dUfdt_func(t);

    %tack the two ends onto either side so every mass can use the same
    %u_{i-1} - 2u_i + u_{i+1} formula. the wall end is always 0 and the
    %right end follows the actuator, which takes care of the boundary
    %conditions (eqns 11 and 12) without any special cases
    U_padded = [0;U;Uf];
    dUdt_padded = [0;dUdt;dUfdt];

    %discrete laplacian for all n masses at once, entry i is
    %(left neighbor) - 2*(mass i) + (right neighbor)
    lap_U = U_padded(1:n) - 2*U_padded(2:(n+1)) + U_padded(3:(n+2));
    lap_dUdt = dUdt_padded(1:n) - 2*dUdt_padded(2:(n+1)) + dUdt_padded(3:(n+2));

    %compute acceleration
    %each mass is M/n, so divide both sides of eqn 10 by M/n
    d2Udt2 = (n/M)*((Tf/dx)*lap_U + (c/dx)*lap_dUdt);

    %assemble state derivative
    dVdt = [dUdt;d2Udt2];
end
