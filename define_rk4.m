%Butcher tableau for the classic 4th order RK method (RK4)
%4 stages, local truncation error scales like h^5
%INPUTS:
%   None
%OUTPUTS:
%RK4: struct that contains the Butcher tableau
%   RK4.A: matrix of a_{ij} values
%   RK4.B: vector of b_i values
%   RK4.C: vector of c_i values
function RK4 = define_rk4()
    RK4 = struct();
    RK4.C = [0, 1/2, 1/2, 1];
    RK4.B = [1/6, 1/3, 1/3, 1/6];
    RK4.A = [0,   0,   0, 0;...
             1/2, 0,   0, 0;...
             0,   1/2, 0, 0;...
             0,   0,   1, 0];
end
