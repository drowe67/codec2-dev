% sd.m
% David Rowe Aug 2012
% Plots the spectal distorion between twofiles of LPCs.  Used for LSP 
% quantisation tuning.

function sd(raw_filename, dump_file_prefix, f)

  graphics_toolkit ("gnuplot");

  e_filename = sprintf("%s_E.txt", dump_file_prefix);
  e = load(e_filename);
  ak1_filename = sprintf("%s_ak.txt", dump_file_prefix);
  ak2_filename = sprintf("%s_ak_.txt", dump_file_prefix);
  ak1 = load(ak1_filename);
  ak2 = load(ak2_filename);

  [ak1_r, ak1_c] = size(ak1);
  [ak2_r, ak2_c] = size(ak1);

  frames = max([ak1_r ak2_r]); printf("%d frames\n", frames);
  sd = zeros(1,frames);
  Ndft = 512;
  A1 = zeros(frames, Ndft);
  A2 = zeros(frames, Ndft);

  % initial helicopter view of all frames

  spec_err = zeros(1, Ndft);
  for i = 1:frames
    A1(i,:) = -20*log10(abs(fft(ak1(i,:),Ndft)));
    A2(i,:) = -20*log10(abs(fft(ak2(i,:),Ndft)));
    sd(i) = sum((A1(i,:) - A2(i,:)).^2)/Ndft;
    spec_err += (A1(i,:) - A2(i,:)).^2;
  end
  spec_err /= frames;
  printf("sd av %3.2f dB*dB\n", sum(sd)/frames);

  % work out worst frames with sig energy
  
  ind = find(e < 0);
  sd(ind) = 0;
  [largest largest_ind] = sort(sd,"descend");
  printf("largest SD frames....: %3.2f\n", largest(1:5));
  printf("largest SD frames ind: %d\n", largest_ind(1:5));
  
  figure(1);
  clf;
  subplot(211)
  fs=fopen(raw_filename,"rb");
  s = fread(fs,Inf,"short");
  plot(s);
  axis([1 length(s) -20E3 20E3])
  subplot(212)
  [a b c] = plotyy(1:frames, sd, 1:frames, e);
  axis(a, [1 frames 0 10])

  lsp1_filename = sprintf("%s_lsp.txt", dump_file_prefix);
  lsp2_filename = sprintf("%s_lsp_.txt", dump_file_prefix);
  lsp1 = load(lsp1_filename);
  lsp2 = load(lsp2_filename);

  mel_filename = sprintf("%s_mel.txt", dump_file_prefix);
  mel = load(mel_filename);

  weights_filename = sprintf("%s_weights.txt", dump_file_prefix); 
  if file_in_path(".",weights_filename)
    weights = load(weights_filename);
  end

  figure(4)
  plot(e,sd,'+')
  axis([0 50 0 10])
  xlabel('LPC energy dB')
  ylabel('SD dB*dB');

  figure(5)
  subplot(211)
  ind = find(e > 0);
  hist(sd(ind))
  subplot(212)
  plot((1:Ndft)*8000/Ndft, spec_err)
  axis([300 3000 0 max(spec_err)])

  % now enter single step mode so we can analyse each frame
  k = ' ';
  largest_mode = 0;
  do 
    if largest_mode
      fr = largest_ind(f);
    else
      fr = f; 
    endif

    figure(2);
    clf;
    plot((4000/pi)*lsp1((fr-2:fr+2),:));
    hold on;
    plot((4000/pi)*lsp2((fr-2:fr+2),:),'+-');
    hold off;

    figure(3);
    clf;

    plot((1:Ndft/2)*4000/(Ndft/2), A1(fr,1:(Ndft/2)),";A1;r");
    axis([1 4000 -20 40]);
    hold on;
    plot((1:Ndft/2)*4000/(Ndft/2), A2(fr,1:(Ndft/2)),";A2;");
    if file_in_path(".",weights_filename)
      plot(lsp1(fr,:)*4000/pi, weights(fr,:),";weights;g+");
    end

    printf("\n");
    for l=1:6
        plot([lsp1(fr,l)*4000/pi lsp1(fr,l)*4000/pi], [0  -10], 'r');
        plot([lsp2(fr,l)*4000/pi lsp2(fr,l)*4000/pi], [-10 -20], 'b');
        plot([mel(fr,l) mel(fr,l)], [0 10], 'g');
        printf("%d  ", mel(fr,l));
    endfor
    printf("\n");
    
    plot(0,0,';lsp1;r');
    plot(0,0,';lsp2;b');
    sd_str = sprintf(";sd %3.2f dB*dB;", sd(fr));
    plot(0,0,sd_str);
   
    hold off;

    % interactive menu

    printf("\rframe: %d  menu: n-next  b-back  l-largest mode  q-quit", fr);
    fflush(stdout);
    k = kbhit();
    if (k == 'n')
      f = f + 1;
    endif
    if (k == 'b')
      f = f - 1;
    endif

    if (k == 'l')
      if largest_mode
        largest_mode = 0;
      else
        largest_mode = 1;
        f = 1;
      endif
    endif

  until (k == 'q')
  printf("\n");

endfunction

