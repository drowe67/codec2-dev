% dec_pager.m
%
% Interactive Octave script to explore variable decimation points for low rate codecs
%
% Usage:
%   Make sure codec2-dev is compiled with the -DDUMP option - see README for
%    instructions.
%   ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --dump hts1a
%   $ cd ~/codec2-dev/octave
%   octave:14> dec_pager("../raw/hts1a.raw", "../build_linux/src/hts1a")

function dec_pager(rawfilename, samname, f=1, varargin)
  more off;

  Fs = 8000; rate_K_sample_freqs_kHz = [0.1:0.1:4]; K = length(rate_K_sample_freqs_kHz);
  Nsam = 80;
   
  % Number of frames (time block) we consider for sample points.  There are typically
  % Nt/dec = 4 samples points per block
  
  Nt = 16;

  % read in optional vector that defines sampling points for
  % decimation, each entry is the frame number of a sample e.g.
  % [1 8 10 12 ....
  
  decsamfile_en = 0;
  ind = arg_exists(varargin, "decsamfile");
  if ind
    dec_filename = varargin{ind+1};
    decvec = load(dec_filename); Ndv = length(decvec);
    decsamfile_en = 1;
  end
  
  % load up raw samples and text files dumped from c2sim -----------------------

  fs = fopen(rawfilename,"rb");
  Sraw = fread(fs,Inf,"short");

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
    % compute rate K surface and mean energy

    rate_K_surf = resample_const_rate_f(model(f:f+Nt-1,:), rate_K_sample_freqs_kHz, Fs);
    rate_k_mean = mean(rate_K_surf,2);

    offset = 2;
    st = (f-1-offset)*Nsam + 1; en = st + Nt*Nsam - 1;
    s = Sraw(st:en);
    figure(1); clf;
    subplot(211); plot(st:en, s); axis([st en -20000 20000]);
    subplot(212); plot(f:(f+Nt-1), rate_k_mean); axis([f f+Nt-1 0 80]);
    if decsamfile_en
      hold on; stem(decvec, 60*ones(1,Ndv)); hold off;
    end
    
    % plots ----------------------------------
  
    figure(2); clf;
    mesh(rate_K_surf);
    axis([1 K 1 Nt 0 80])
    
    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back", f);
    fflush(stdout);
    k = kbhit();

    if (k == 'n') && (f+2*Nt <= frames)
      f = f + Nt;
    endif
    if k == 'b'
      f = f - Nt;
    endif
  until (k == 'q')
  printf("\n");

endfunction

 
function ind = arg_exists(v, str) 
   ind = 0;
   for i=1:length(v)
      if strcmp(v{i}, str)
        ind = i;
      end     
    end
endfunction


