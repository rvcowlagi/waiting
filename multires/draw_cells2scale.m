function draw_cells2scale(V, myAxes, color_border, width_border, ...
	txt_cell, cells2draw, color_face, nomcellsize, wksp)


if numel(myAxes)
	axes(myAxes); hold on;
end
if (nargin <= 5) || (numel(cells2draw) == 0)
	cells2draw = 1:size(V, 1);
end

if width_border > 0
	for m = cells2draw
		if (nargin <= 7)
			rect_posn = [V(m,1) V(m,2) V(m,3) V(m,3)];
		else
			rect_posn = [ (-wksp + nomcellsize*V(m,1:2)) nomcellsize*[V(m,3) V(m,3)] ];
		end

		if (nargin > 6) && ((numel(color_face) ~= 1) || (color_face ~= 0))
			rectangle('Position', rect_posn, ...
				'EdgeColor', color_border, 'LineWidth', width_border, ...
				'FaceColor', color_face);
		else
			rectangle('Position', rect_posn, ...
				'EdgeColor', color_border, 'LineWidth', width_border);
		end
		if txt_cell
			text( -wksp + nomcellsize*(V(m,1) + 0.1*V(m,3)), ...
				-wksp + nomcellsize*(V(m,2) + 0.15*V(m,3)), num2str(m), ...
				'FontName', 'Consolas', 'FontSize', txt_cell, ...
				'FontWeight', 'bold', 'Color', 'r');
        end
	end
end
drawnow;