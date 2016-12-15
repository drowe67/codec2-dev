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
  quant_en = 0; pf_en = 0;
  melvq;

  K=20; load train_120_vq; m=5; 

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
    plot(s);
    axis([1 length(s) -20000 20000]);
    title('Time Domain Speech');

    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    % plots for mel sampling

    [rate_K_vec rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model(f,:), K);

    mean_f = mean(rate_K_vec);
    rate_K_vec_no_mean = rate_K_vec - mean_f;
    
    if quant_en == 2
      [res rate_K_vec_no_mean_ ind] = mbest(train_120_vq, rate_K_vec_no_mean, m);
    else
      rate_K_vec_no_mean_ = rate_K_vec_no_mean;
    end

    if pf_en
      % pf, needs some energy equalisation, does gd things for hts1a
      rate_K_surface_ *= 1.2;
    end

    rate_K_vec_ = rate_K_vec_no_mean_ + mean_f;
    [model_ AmdB_] = resample_rate_L(model(f,:), rate_K_vec_, rate_K_sample_freqs_kHz);

    % plots ----------------------------------

    figure(2);
    clf;
    title('Frequency Domain 1');

    axis([1 4000 -20 80]);
    hold on;
    plot((1:L)*Wo*4000/pi, AmdB,";Am;b+-");
    plot(rate_K_sample_freqs_kHz*1000, rate_K_vec, ';rate K mel;g+-');
    if quant_en >= 1
      plot((1:L)*Wo*4000/pi, AmdB_,";Am quant;k+-");
    end
    if quant_en == 2
      plot(rate_K_sample_freqs_kHz*1000, rate_K_vec_, ';rate K mel quant;r+-');   
    end

    hold off;

    figure(3);
    clf;
    title('Frequency Domain 2');
    axis([1 4000 -80 80]);
    hold on;
    plot((1:L)*Wo*4000/pi, AmdB,";Am;b+-");
    plot(rate_K_sample_freqs_kHz*1000, rate_K_vec_no_mean, ';rate K mel no mean;g+-');
    if quant_en == 2
    plot(rate_K_sample_freqs_kHz*1000, rate_K_vec_no_mean_, ';rate K mel no mean quant;r+-');
    end
    hold off;

    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back  q-quit  m-quant_en[%d]", f, quant_en);
    fflush(stdout);
    k = kbhit();

    if (k == 'm')
      quant_en++;
      if quant_en > 2
        quant_en = 0;
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

 
#{ Piecewise model stuff, organise later if rqd
    [maskdB Am_freqs_kHz] = mask_model(AmdB, Wo, L);
    AmdB_ = maskdB;
    [mx mx_ind] = max(AmdB_);
    AmdB_(mx_ind) += 6;
   

    if quant_en
      [AmdB_ residual fvec fvec_] = piecewise_model(AmdB, Wo, vq, 1);
    else
      [AmdB_ residual fvec] = piecewise_model(AmdB, Wo);
    end
    fvec
#}
