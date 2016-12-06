% newamp1_fbf.m
%
% Copyright David Rowe 2016
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Interactive Octave script to explore frame by frame operation of new amplitude
% modelling model.
%
% Usage:
%   Make sure codec2-dev is compiled with the -DDUMP option - see README for
%    instructions.
%   ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --dump hts1a
%   $ cd ~/codec2-dev/octave
%   octave:14> newamp1_fbf("../build_linux/src/hts1a",50)



function newamp1_fbf(samname, f=10)
  newamp;
  more off;
  quant_en = 0;

  load vq;

  % load up text files dumped from c2sim ---------------------------------------

  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);
  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);

  % Keyboard loop --------------------------------------------------------------

  k = ' ';
  do 
    figure(1);
    clf;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    size(s);
    plot(s);
    axis([1 length(s) -20000 20000]);
    title('Time Domain Speech');

    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    #{   
    [maskdB Am_freqs_kHz] = mask_model(AmdB, Wo, L);
    AmdB_ = maskdB;
    [mx mx_ind] = max(AmdB_);
    AmdB_(mx_ind) += 6;
    #}

    if quant_en
      [AmdB_ residual fvec fvec_] = piecewise_model(AmdB, Wo, vq, 2);
    else
      [AmdB_ residual fvec] = piecewise_model(AmdB, Wo);
    end

    figure(2);
    clf;
    title('Frequency Domain');

    axis([1 4000 -20 80]);
    hold on;
    plot((1:L)*Wo*4000/pi, AmdB,";Am;r+-");
    plot(Am_freqs_kHz*1000, AmdB_, ';model;c');
    plot(fvec*1000, 60*ones(1,4), ';fvec;go');
    if quant_en
      plot(fvec_*1000, 60*ones(1,4), ';fvec q;ro');
    end

    hold off;

    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back  q-quit  m-quant_en", f);
    fflush(stdout);
    k = kbhit();

    if (k == 'm')
      if quant_en
        quant_en = 0;
      else
        quant_en = 1; 
      end
    endif
    if (k == 'n')
      f = f + 1;
    endif
    if (k == 'b')
      f = f - 1;
    endif
  until (k == 'q')
  printf("\n");

endfunction

 
