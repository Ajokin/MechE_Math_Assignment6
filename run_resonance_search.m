close all;
% for n = 2 and n = 3 masses, find the actuator frequencies that excite
% each mode. modal analysis says where the resonances should be, then a
% sweep of omega through the actual simulation checks that the big
% responses land right on them

% figures for the report, saved one folder up so they don't end up in the repo
fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Assignment 6 Report Figures');
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% first make sure the matrix version of the motion equation (eqn 22)
% gives the same accelerations as the rate function at some random state,
% otherwise the modes won't mean anything for the simulation
test_params = define_string_params(4,0.1,2);
V_test = rand(8,1);
t_test = 0.7;
dVdt_test = string_rate_func01(t_test,V_test,test_params);
[M_mat,K_mat,C_mat,B1,B2] = construct_2nd_order_matrices(test_params);
d2Udt2_matrix = M_mat\(-C_mat*V_test(5:8) - K_mat*V_test(1:4) ...
    + B1*test_params.Uf_func(t_test) + B2*test_params.dUfdt_func(t_test));
fprintf('max difference, rate function vs matrix form = %.3e\n',max(abs(dVdt_test(5:8)-d2Udt2_matrix)));

RK4 = define_rk4();
h_ref = 0.02;
tspan = [0,40];
amplitude_Uf = 0.01;

for num_masses = [2,3]
    string_params = define_string_params(num_masses,0,0);
    [Ur_mat,omega_n] = compute_string_modes(string_params);

    % closed form for n equal masses on a string, from the discrete mode
    % shapes link in the pdf, just to double check eig
    k = string_params.Tf/string_params.dx;
    m = string_params.M/num_masses;
    omega_closed_form = 2*sqrt(k/m)*sin((1:num_masses)'*pi/(2*(num_masses+1)));

    for i = 1:num_masses
        fprintf('n = %d, mode %d: omega = %.4f (closed form %.4f), Ur = [%s]\n', ...
            num_masses,i,omega_n(i),omega_closed_form(i),strtrim(sprintf('%.3f ',Ur_mat(:,i))));
    end

    % sweep from well under the first mode to past the last one, and add
    % the exact resonant frequencies in so the peaks don't get skipped
    omega_list = sort([linspace(0.5*omega_n(1),1.25*omega_n(end),150),omega_n']);
    peak_ratio_list = zeros(size(omega_list));

    for j = 1:length(omega_list)
        string_params = define_string_params(num_masses,amplitude_Uf,omega_list(j));
        my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);
        V0 = zeros(2*num_masses,1);
        [tlist,Vlist] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0,h_ref,RK4);

        % biggest swing of any mass over the whole run, compared to how
        % far the actuator moves
        peak_ratio_list(j) = max(max(abs(Vlist(:,1:num_masses))))/amplitude_Uf;
    end

    for i = 1:num_masses
        fprintf('n = %d: at omega_%d, peak |u| / A = %.2f\n',num_masses,i,peak_ratio_list(omega_list == omega_n(i)));
    end

    fig = figure('Color','w','Position',[100,100,1100,650]);
    hold on; box on; grid on;
    set(gca,'YScale','log');

    sweep_plot = plot(omega_list,peak_ratio_list,'-','Color',[0 .25 .85],'LineWidth',2.5);
    for i = 1:num_masses
        mode_line = xline(omega_n(i),'--',sprintf('$\\omega_%d$ = %.3f',i,omega_n(i)), ...
            'Color',[.9 .15 .05],'LineWidth',2,'Interpreter','latex','FontSize',16, ...
            'LabelOrientation','horizontal','LabelVerticalAlignment','bottom');
    end
    % ratio of 1 means the masses move as much as the actuator does
    one_line = yline(1,'k:','LineWidth',1.5);

    xlim([omega_list(1),omega_list(end)]);
    % extra room on top (log scale) so the legend sits above the peaks
    % instead of covering the first one
    ylim([0.5*min(peak_ratio_list),5*max(peak_ratio_list)]);
    xlabel('Actuator Frequency, $\omega$ (-)','Interpreter','latex');
    ylabel('Peak Displacement Ratio, $\max|u_i|/A$ (-)','Interpreter','latex');
    title(sprintf('Frequency Sweep of the Simulation, %d Masses ($A$ = %g, Starting From Rest)',num_masses,amplitude_Uf), ...
        'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex','FontSize',16);
    legend([sweep_plot,mode_line,one_line],{'Simulation','Resonant frequencies from modal analysis, $\omega_r$', ...
        'String moves as much as the actuator'},'Interpreter','latex','Location','northwest','FontSize',13);

    exportgraphics(fig,fullfile(fig_dir,sprintf('frequency_sweep_%d_masses.png',num_masses)),'Resolution',300);
end
