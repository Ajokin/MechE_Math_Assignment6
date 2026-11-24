%for three of the mode shapes, drive the actuator right at that mode's
%resonant frequency (starting from rest) and compare the simulation to
%the predicted mode shape. one video per mode, top subplot is the
%simulation and the bottom one is the mode shape from the eigenvector
%INPUTS:
%mode_list: (optional) which modes to make videos of, defaults to 1:3.
%   each video takes a while to render so this lets them run separately
function run_mode_shape_animation(mode_list)
    if nargin < 1
        mode_list = 1:3;
    end
    close all;

    % videos go one folder up with the other report figures
    fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Assignment 6 Report Figures');
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end

    % same string as the basic animation
    num_masses = 5;
    string_params = define_string_params(num_masses,0,0);

    % mode shapes are the eigenvectors of the generalized eigenvalue
    % problem K*Ur = omega^2*M*Ur (eqn 26)
    [Ur_mat,omega_n] = compute_string_modes(string_params);
    for i = 1:num_masses
        fprintf('mode %d: omega = %.4f, Ur = [%s]\n',i,omega_n(i),strtrim(sprintf('%.3f ',Ur_mat(:,i))));
    end

    %list of x points (including the two endpoints)
    xlist = linspace(0,string_params.L,num_masses+2);

    RK4 = define_rk4();
    frame_rate = 60;
    h_ref = 1/(4*frame_rate);
    tspan = [0,30];
    t_frames = (tspan(1):1/frame_rate:tspan(2))';

    % picked so all three videos end up swinging about the same amount (peak
    % around 0.15). the higher modes have more damping, zeta = c*omega/(2*Tf),
    % so they don't build up as much and need a bigger push
    amplitude_list = [0.006, 0.009, 0.012];

    for mode_num = mode_list
        Ur = Ur_mat(:,mode_num);
        omega_r = omega_n(mode_num);
        amplitude_Uf = amplitude_list(mode_num);

        string_params = define_string_params(num_masses,amplitude_Uf,omega_r);
        my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);

        %string starts at rest, all the masses at zero and not moving
        V0 = zeros(2*num_masses,1);
        [tlist,Vlist] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0,h_ref,RK4);
        U_list = Vlist(:,1:num_masses);
        U_peak = max(max(abs(U_list)));

        % check how much of the simulation is actually in this mode. all the
        % masses are the same so the modes are orthogonal, and the amount of
        % mode r in U is just a dot product, q(t) = Ur'*U/(Ur'*Ur)
        q_list = U_list*Ur/(Ur'*Ur);
        leftover = U_list - q_list*Ur';
        fprintf('mode %d: peak |u| / A = %.2f, biggest part not in the mode shape / peak = %.3f\n', ...
            mode_num,U_peak/amplitude_Uf,max(max(abs(leftover)))/U_peak);

        % predicted motion straight from the eigenvector. driving right at
        % resonance the response lags the actuator by 90 degrees, so it goes
        % like sin(omega_r*t) instead of cos. the actuator pushes on the last
        % mass, so the sign follows Ur(n). the eigenvector's size doesn't
        % mean anything, so it's scaled to the biggest swing of the
        % simulation just so the two plots are the same size
        mode_sign = sign(Ur(end));
        U_pred_amp = U_peak*mode_sign*Ur;

        % evenly spaced frame times so the video plays at a constant speed
        V_frames = interp1(tlist,Vlist,t_frames);

        % both subplots get the same fixed axes so they're easy to compare
        % (extra room so the time stamp fits under the masses)
        y_max = 1.6*U_peak;

        % dashed envelope of the predicted mode shape, drawn on both
        % subplots. the simulated string should grow into it
        envelope = U_peak*[0;Ur;0];

        input_fname = fullfile(fig_dir,sprintf('mode_%d_string_animation.mp4',mode_num));
        writerObj = VideoWriter(input_fname,'MPEG-4');
        writerObj.FrameRate = frame_rate;
        writerObj.Quality = 100;
        open(writerObj);

        %initialize the current figure and save as object
        fig1 = figure(1);
        %this is the critical line that forces the video to be high quality
        set(fig1,'units','pixels','position',[0 0 1920 1080]);
        set(fig1,'color','w');
        %same as string_simulation, print at the paper size so the frames are
        %exactly 1920 x 1080 even on the small laptop screen
        set(fig1,'PaperUnits','inches','PaperPosition',[0 0 19.2 10.8]);
        tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
        title(tl,sprintf('Mode %d: Simulation vs. Prediction ($\\omega = \\omega_%d$ = %.3f, $A$ = %g)', ...
            mode_num,mode_num,omega_r,amplitude_Uf),'Interpreter','latex','FontSize',28);

        %top: the actual simulation, end driven by the actuator
        nexttile;
        sim_drawing = initialize_string_drawing(string_params,y_max);
        env_plot = plot(xlist,envelope,'--','Color',[.9 .15 .05],'LineWidth',2);
        plot(xlist,-envelope,'--','Color',[.9 .15 .05],'LineWidth',2);
        %full label is too tall for half the figure
        ylabel('Displacement, $u$ (-)','Interpreter','latex');
        title(sprintf('Simulation (%d Masses), Starting From Rest',num_masses), ...
            'Interpreter','latex','FontSize',24);

        %bottom: predicted mode shape from the eigenvector, both ends held
        %at 0 like in the eigenvalue problem
        nexttile;
        mode_drawing = initialize_string_drawing(string_params,y_max);
        plot(xlist,envelope,'--','Color',[.9 .15 .05],'LineWidth',2);
        plot(xlist,-envelope,'--','Color',[.9 .15 .05],'LineWidth',2);
        ylabel('Displacement, $u$ (-)','Interpreter','latex');
        title(sprintf('Predicted Mode %d Shape From the Eigenvector, $U_%d\\sin(\\omega_%d t)$', ...
            mode_num,mode_num,mode_num),'Interpreter','latex','FontSize',24);

        %one legend for both plots, off to the right so it can't cover anything
        lgd = legend([sim_drawing.string_line,sim_drawing.masses,env_plot,sim_drawing.ends,sim_drawing.actuator], ...
            {'String','Masses','Predicted mode shape envelope','String endpoints, $u_0$ and $u_f(t)$', ...
            'Linear actuator, $u_f(t)$'},'Interpreter','latex','FontSize',18);
        lgd.Layout.Tile = 'east';

        for i = 1:length(t_frames)
            update_string_drawing(sim_drawing,V_frames(i,1:num_masses)', ...
                string_params.Uf_func(t_frames(i)),t_frames(i));
            update_string_drawing(mode_drawing,U_pred_amp*sin(omega_r*t_frames(i)),0,t_frames(i));
            drawnow;

            %capture a frame (what is currently plotted)
            current_frame = print(fig1,'-RGBImage','-r100','-opengl');

            %write the frame to the video
            writeVideo(writerObj,current_frame);
        end

        %must call close after all frames are written to save the video
        close(writerObj);
        close(fig1);
    end
end
