% fdmdv_demod_c.m
%
% Plots Octave dump file information from C FDMDV demodulator program,
% to give a similar set of plots to fdmdv_demod.m.  Useful for off
% line analysis of demod performance.
%
% Copyright David Rowe 2012
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%

function fdmdv_demod_c(dumpfilename, bits)

  fdmdv; % include modem code
  frames = bits/(Nc*Nb);

  load(dumpfilename);

  % ---------------------------------------------------------------------
  % Plots
  % ---------------------------------------------------------------------

  xt = (1:frames)/Rs;
  secs = frames/Rs;

  figure(1)
  clf;
  plot(real(rx_symbols_log_c(1:Nc+1,15:frames)),imag(rx_symbols_log_c(1:Nc+1,15:frames)),'+')
  axis([-2 2 -2 2]);
  title('Scatter Diagram');

  figure(2)
  clf;
  subplot(211)
  plot(xt, rx_timing_log_c(1:frames))
  title('timing offset (samples)');
  subplot(212)
  plot(xt, foff_log_c(1:frames))
  hold on;
  plot(xt, coarse_fine_log_c(1:frames)*75, 'r');
  hold off;
  title('Freq offset (Hz)');
  grid

endfunction
