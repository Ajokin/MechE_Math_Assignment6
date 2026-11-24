%Updates the plot objects that show the vibrating string
%for the current displacement of the masses and the actuator
%INPUTS:
%string_drawing: a struct containing all the plotting objects
%       (made by initialize_string_drawing)
%U: n x 1 column of the transverse displacements of the masses
%Uf: current displacement of the actuator end, uf(t)
%t: the current time (for the time stamp)
function update_string_drawing(string_drawing,U,Uf,t)
    %tack the two ends onto the masses, the wall end is always at 0
    %and the right end follows the actuator
    U_plot = [0;U;Uf];

    set(string_drawing.string_line,'ydata',U_plot);
    set(string_drawing.masses,'ydata',U);
    set(string_drawing.ends,'ydata',[0,Uf]);
    set(string_drawing.rod,'ydata',[string_drawing.rod_bottom,Uf]);

    set(string_drawing.time_text,'string',sprintf('Time, $t$ = %.2f (-)',t));
end
