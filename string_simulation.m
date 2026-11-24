%runs the vibrating string simulation (masses on a taut string, driven by
%the linear actuator on the right end) and saves the basic animation
function string_simulation()
    close all;

    %works for any number of masses, 5 looks nice
    num_masses = 5;

    %no forcing yet, just need the resonant frequencies of this string
    string_params = define_string_params(num_masses,0,0);
    [~,omega_n] = compute_string_modes(string_params);

    %drive it right at the 2nd resonant frequency so the masses build up
    %into the mode 2 shape. the actuator only moves a tiny bit so you can
    %see the string swinging way more than the end is being pushed
    amplitude_Uf = 0.009;
    omega_Uf = omega_n(2);
    string_params = define_string_params(num_masses,amplitude_Uf,omega_Uf);

    %load string_params into rate function
    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);

    %initial conditions, string starts flat and not moving
    U0 = zeros(num_masses,1);
    dUdt0 = zeros(num_masses,1);
    V0 = [U0;dUdt0];

    tspan = [0,30];

    %run the integration
    %fixed step RK4, 4 steps per video frame. way under the stability
    %limit (h*omega_max is about 0.04 for 5 masses)
    frame_rate = 60;
    RK4 = define_rk4();
    h_ref = 1/(4*frame_rate);
    [tlist,Vlist,h_avg,num_evals] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0,h_ref,RK4);
    fprintf('%d steps, h_avg = %.5f, %d rate function calls\n',length(tlist)-1,h_avg,num_evals);

    %how much bigger the string swing gets than the actuator motion
    U_list = Vlist(:,1:num_masses);
    U_peak = max(max(abs(U_list)));
    fprintf('omega = %.4f, peak |u| = %.4f, peak |u| / A = %.2f\n',omega_Uf,U_peak,U_peak/amplitude_Uf);

    %1 time unit per second of video so it plays back in real time
    t_frames = (tspan(1):1/frame_rate:tspan(2))';
    V_frames = interp1(tlist,Vlist,t_frames);

    %axes are set from the biggest swing in the whole run so they stay
    %put for the entire video, 40% extra so nothing touches the edge
    y_max = 1.4*U_peak;

    %saves the video one folder up (right outside the repo) with the
    %other report figures
    fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Assignment 6 Report Figures');
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    input_fname = fullfile(fig_dir,'vibrating_string_basic_animation.mp4');

    %mp4 instead of avi so the file is small enough to upload quick
    writerObj = VideoWriter(input_fname,'MPEG-4');
    writerObj.FrameRate = frame_rate;
    writerObj.Quality = 100;
    open(writerObj);

    %initialize the current figure and save as object
    fig1 = figure(1);

    %this is the critical line that forces the video to be high quality
    %1920 x 1080 so it uploads to youtube as 1080p
    set(fig1,'units','pixels','position',[0 0 1920 1080]);
    set(fig1,'color','w');

    %the laptop screen is only 1536 x 864 (125% scaling) so windows
    %squishes the figure and getframe came out 2400 x 990. print renders
    %the figure off screen at the paper size instead, 19.2 x 10.8 inches
    %at 100 dpi is exactly 1920 x 1080 on any computer. it's a lot slower
    %than getframe (about 1 sec a frame) so the video takes a while
    set(fig1,'PaperUnits','inches','PaperPosition',[0 0 19.2 10.8]);

    string_drawing = initialize_string_drawing(string_params,y_max);

    title({sprintf('Vibrating String Animation (%d Masses)',num_masses), ...
        sprintf('Actuator Driven at Mode 2: $u_f(t) = A\\cos(\\omega t)$, $A$ = %g, $\\omega = \\omega_2$ = %.3f', ...
        amplitude_Uf,omega_Uf)},'Interpreter','latex','FontSize',28);

    legend([string_drawing.string_line,string_drawing.masses,string_drawing.ends, ...
        string_drawing.wall,string_drawing.actuator,string_drawing.equilibrium], ...
        {'String','Masses, $m_i = M_{total}/n$','String endpoints, $[0,0]$ and $[L,u_f(t)]$', ...
        'Fixed wall','Linear actuator, $u_f(t)$','Equilibrium, $u = 0$'}, ...
        'Interpreter','latex','Location','southoutside','NumColumns',3,'FontSize',18);

    for i = 1:length(t_frames)
        U = V_frames(i,1:num_masses)';
        Uf = string_params.Uf_func(t_frames(i));
        update_string_drawing(string_drawing,U,Uf,t_frames(i));
        drawnow;

        %capture a frame (what is currently plotted)
        %opengl is about twice as fast as the default for print
        current_frame = print(fig1,'-RGBImage','-r100','-opengl');

        %write the frame to the video
        writeVideo(writerObj,current_frame);
    end

    %must call close after all frames are written to save the video
    close(writerObj);
end
