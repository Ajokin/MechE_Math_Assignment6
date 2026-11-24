%Defines the parameters for the vibrating string. everything except the
%number of masses and the actuator motion stays the same for every run, so
%the simulation and the modal analysis scripts all use the exact same string
%INPUTS:
%num_masses: number of masses on the string
%amplitude_Uf: amplitude of the actuator motion, uf(t) = A*cos(omega*t)
%omega_Uf: frequency of the actuator motion
%OUTPUTS:
%string_params: a struct containing the system parameters describing the string
%   string_params.n: number of masses
%   string_params.M: total mass attached to the string
%   string_params.Uf_func: function describing motion of end point
%   string_params.dUfdt_func: time derivative of Uf
%   string_params.Tf: %tension in string
%   string_params.L: %length of string
%   string_params.c: %damping coefficient
%   string_params.dx: %horizontal spacing between masses
%   string_params.amplitude_Uf: A (kept so the titles can print it)
%   string_params.omega_Uf: omega (kept so the titles can print it)
function string_params = define_string_params(num_masses,amplitude_Uf,omega_Uf)
    %total mass stays fixed, so more masses means each one is lighter
    total_mass = 1;
    tension_force = 1;
    string_length = 1;

    %small number
    %just enough damping to keep the resonant runs from blowing up (eqn 10)
    %tried 0.005 first but mode 3 only got to about 4x the actuator
    %amplitude, at 0.001 all three modes get past 12x in 30 time units
    damping_coeff = 0.001;

    %n masses split the string into n+1 even pieces
    dx = string_length/(num_masses+1);

    Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
    dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);

    %generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    string_params.amplitude_Uf = amplitude_Uf;
    string_params.omega_Uf = omega_Uf;
end
