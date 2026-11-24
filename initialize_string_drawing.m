%Creates the plotting objects for the vibrating string animation
%(the string, the masses, the wall on the left and the linear actuator on
%the right) in whatever axes are current, so the same drawing can go in a
%full figure or in one of the subplots
%INPUTS:
%string_params: a struct containing the system parameters describing the string
%y_max: biggest displacement the axes need to show, the y axis is fixed
%       at [-y_max, y_max] for the whole video
%OUTPUTS:
%string_drawing: a struct containing all the plotting objects
%       string_drawing.string_line is the line through every point on the string
%       string_drawing.masses is the plot of the n masses
%       string_drawing.ends is the two endpoints, [0,0] and [L,uf(t)]
%       string_drawing.wall is the patch for the wall on the left
%       string_drawing.actuator is the patch for the actuator body
%       string_drawing.rod is the actuator rod that moves with uf(t)
%       string_drawing.rod_bottom is where the rod comes out of the actuator
%       string_drawing.equilibrium is the dashed line at u = 0
%       string_drawing.time_text shows the current time
function string_drawing = initialize_string_drawing(string_params,y_max)
    L = string_params.L;
    n = string_params.n;

    %list of x points (including the two endpoints)
    xlist = linspace(0,L,n+2);

    hold on; box on; grid on;
    axis equal;
    %fixed limits so the axes don't shift around during the video
    %a little room on the left for the wall and on the right for the actuator
    axis([-0.08*L, 1.08*L, -y_max, y_max]);

    string_drawing = struct();

    %background stuff goes first so the string and masses draw on top
    string_drawing.equilibrium = plot([0,L],[0,0],'--','Color',[.55 .55 .55],'LineWidth',1.5);

    %wall on the left, the string is pinned to its right face
    wall_w = 0.04*L;
    wall_h = 0.5*y_max;
    string_drawing.wall = patch([-wall_w,0,0,-wall_w],[-wall_h,-wall_h,wall_h,wall_h], ...
        [.35 .35 .35],'EdgeColor','k','LineWidth',2);

    %actuator body sits below the string on the right, the rod comes out
    %the top of it and holds the end of the string at uf(t)
    act_w = 0.03*L;
    act_top = -0.45*y_max;
    act_bottom = -0.9*y_max;
    string_drawing.actuator = patch([L-act_w/2,L+act_w/2,L+act_w/2,L-act_w/2], ...
        [act_bottom,act_bottom,act_top,act_top],[.8 .8 .8],'EdgeColor','k','LineWidth',2);
    string_drawing.rod_bottom = act_top;
    string_drawing.rod = plot([L,L],[act_top,0],'-','Color',[.45 .45 .45],'LineWidth',6);

    string_drawing.string_line = plot(xlist,zeros(1,n+2),'k-','LineWidth',2.5);

    %gray masses like the diagram in the pdf
    string_drawing.masses = plot(xlist(2:(n+1)),zeros(1,n),'o','MarkerSize',24, ...
        'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','k','LineWidth',2);

    string_drawing.ends = plot([0,L],[0,0],'o','MarkerSize',11, ...
        'MarkerFaceColor',[.9 .15 .05],'MarkerEdgeColor','k','LineWidth',1.5);

    %bottom left corner, sitting on the bottom edge so it stays under the
    %masses when they swing down
    string_drawing.time_text = text(0.01*L,-0.97*y_max,'','Interpreter','latex','FontSize',24, ...
        'VerticalAlignment','bottom');

    ax = gca;
    ax.FontSize = 22;
    ax.TickLabelInterpreter = 'latex';
    xlabel('Position Along String, $x$ (-)','Interpreter','latex');
    ylabel('Transverse Displacement, $u$ (-)','Interpreter','latex');
end
