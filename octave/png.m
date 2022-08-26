% Copyright David Rowe 2009
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Replot current plot as a png

function png(pngname)
    print(pngname, "-dpng", "-S800,600");
endfunction
