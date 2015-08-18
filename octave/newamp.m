% Copyright David Rowe 2015
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Octave script to explore new ideas in amplitude modelling.
%
% we don't care about
%  + spectral tilt, in can vary on input and our quantiser shouldnt care.
%    We can vary it on output and the listener won't care 
%  + absolute energy of entire signal
%  + harmonics beneath the masking curve
% we do care about:
%  + clearly defined formant formation

1;

% frame by frame analysis

function newamp_frame(samname, f)
  
  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);

  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);


  Ew_on = 1;
  k = ' ';
  do 
    figure(1);
    clf;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    size(s);
    plot(s);
    axis([1 length(s) -20000 20000]);

    figure(2);
    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);

    % plotting

    plot((1:L)*Wo*4000/pi, AmdB,";Am;r");
    axis([1 4000 0 80]);
    hold on;
    plot((0:255)*4000/256, Sw(f,:),";Sw;");

    [maskdB mask_sample_freqs_kHz] = determine_mask(AmdB, Wo, L);
    non_masked_m = find(maskdB < AmdB);
    plot(mask_sample_freqs_kHz*1000, maskdB, 'g');
    plot(non_masked_m*Wo*4000/pi, AmdB(non_masked_m),'b+');

    hold off;

    % interactive menu

    printf("\rframe: %d  menu: n-next  b-back  p-png  q-quit", f);
    fflush(stdout);
    k = kbhit();
    if (k == 'n')
      f = f + 1;
    endif
    if (k == 'b')
      f = f - 1;
    endif

  until (k == 'q')
  printf("\n");

endfunction


function [maskdB mask_sample_freqs_kHz] = determine_mask(AmdB, Wo, L)

    % calculate and plot masking curve

    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
    maskdB = zeros(1,L);
    for m=1:L
      fmasker_kHz = m*Wo*4/pi;
      maskdB = max(maskdB, schroeder(fmasker_kHz, mask_sample_freqs_kHz) + AmdB(m)); 
    end
end


% process a whole file and write results

function newamp_batch(samname)
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model)
  max_amp = 80;

  Am_out_name = sprintf("%s_am.out", samname);
  fam = fopen(Am_out_name,"wb"); 

  for f=1:frames

    L = min([model(f,2) max_amp-1]);
    Wo = model(f,1);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);
    
    % zero out any harmonics beneath mask

    [maskdB mask_sample_freqs_kHz] = determine_mask(AmdB, Wo, L);
    non_masked_m = find(maskdB < AmdB);
    Am_ = zeros(1,max_amp);
    Am_(2:L) = 10 .^ (maskdB(1:L-1)/20);
    Am_(non_masked_m+1) *=2;
    
    % save to file
    fwrite(fam, Am_, "float32");
  end

  fclose(fam);

endfunction

%
% Masking functions from http://www.perceptualentropy.com/coder.html#C
% Thanks Jon Boley!
%

% Calculate the Schroeder masking spectrum for a given frequency and SPL

function maskdB = schroeder(freq_tone_kHz, mask_sample_freqs_kHz)
  f_kHz = mask_sample_freqs_kHz;
  A = 3.64*(f_kHz).^(-0.8) - 6.5*exp(-0.6*(f_kHz - 3.3).^2) + (10^(-3))*(f_kHz).^4;
  f_Hz = f_kHz*1000;

  % Schroeder Spreading Function
  dz = bark(freq_tone_kHz*1000)-bark(f_Hz);
  maskdB = 15.81 + 7.5*(dz+0.474) - 17.5*sqrt(1 + (dz+0.474).^2);
endfunction


% Converts frequency to bark scale
% Frequency should be specified in Hertz

function b=bark(f)
  b = 13*atan(0.76*f/1000) + 3.5*atan((f/7500).^2); 
endfunction


% plot some masking curves

function plot_masking
  Fs = 8000;
  mask_sample_freqs_kHz = 0.1:0.1:(Fs/1000)/2;
  maskdB = schroeder(1, mask_sample_freqs_kHz)
  figure(1)
  plot(mask_sample_freqs_kHz, maskdB);
endfunction

% Choose one of these to run

more off;
newamp_frame("../build_linux/src/hts2a", 100);
%newamp_batch("../build_linux/src/kristoff");
%plot_masking
