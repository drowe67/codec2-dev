% tfsk.m
% Brady O'Brien 8 January 2016
%
% Octave script to check c port of fsk_horus against the fsk_horus.m
%
% [X] - Functions to wrap around fsk_mod and fsk_demod executables
%     [X] - fsk_mod
%     [X] - fsk_demod
% [X] - Functions to wrap around octave and c implementations, pass
%       same dataset, compare outputs, and give clear go/no-go
%     [X] - fsk_mod_test
%     [X] - fsk_demod_test
% [.] - Port of run_sim and EbNodB curve test battery
% [ ] - Extract and compare more parameters from demod
% [ ] - Run some tests in parallel

fsk_horus_as_a_lib = 1; % make sure calls to test functions at bottom are disabled
fsk_horus_2fsk;  

graphics_toolkit('gnuplot');

fsk_demod_ex_file = "../build/src/fsk_demod";

global mod_pass_fail_maxdiff = 1e-3/50000;

function mod = fsk_mod_c(Fs,Rs,f1,f2,bits)
    %Name of executable containing the modulator
    fsk_mod_ex_file = '../build/src/fsk_mod';
    
    %command to be run by system to launch the modulator
    command = sprintf('%s %d %d %d %d fsk_mod_ut_bitvec fsk_mod_ut_modvec',fsk_mod_ex_file,Fs,Rs,f1,f2);
    
    %save input bits into a file
    bitvecfile = fopen('fsk_mod_ut_bitvec','wb+');
    fwrite(bitvecfile,bits,'uint8');
    fclose(bitvecfile);
    
    %run the modulator
    system(command);
    
    modvecfile = fopen('fsk_mod_ut_modvec','rb');
    mod = fread(modvecfile,'single');
    fclose(modvecfile);
    
endfunction

function bits = fsk_demod_c(Fs,Rs,f1,f2,mod)
    %Name of executable containing the modulator
    fsk_demod_ex_file = '../build/src/fsk_demod';
    modvecfilename = sprintf('fsk_demod_ut_modvec_%d',getpid());
    bitvecfilename = sprintf('fsk_demod_ut_bitvec_%d',getpid());
    
    %command to be run by system to launch the modulator
    command = sprintf('%s %d %d %d %d %s %s',fsk_demod_ex_file,Fs,Rs,f1,f2,modvecfilename,bitvecfilename);
    
    %save modulated input into a file
    modvecfile = fopen(modvecfilename,'wb+');
    fwrite(modvecfile,mod,'single');
    fclose(modvecfile);
    
    %run the modulator
    system(command);
    
    bitvecfile = fopen(bitvecfilename,'rb');
    bits = fread(bitvecfile,'uint8');
    fclose(bitvecfile);
    bits = bits!=0;
    
    %Clean up files
    delete(bitvecfilename);
    delete(modvecfilename);
endfunction

function [dmod,cmod,omod,pass] = fsk_mod_test(Fs,Rs,f1,f2,bits,tname)
    global mod_pass_fail_maxdiff;
    %Run the C modulator
    cmod = fsk_mod_c(Fs,Rs,f1,f2,bits);
    %Set up and run the octave modulator
    states = fsk_horus_init(Fs,Rs);
    states.f1_tx = f1;
    states.f2_tx = f2;
    states.dA = 1;
    states.dF = 0;
    omod = fsk_horus_mod(states,bits);
    
    dmod = cmod-omod;
    (mod_pass_fail_maxdiff*length(dmod))
    max(abs(dmod))
    pass = max(dmod)<(mod_pass_fail_maxdiff*length(dmod))
    if !pass
        printf('Mod failed test %s!\n',tname);
    end
endfunction

function [cbits,obits,pass] = fsk_demod_test(Fs,Rs,f1,f2,mod,tname)
    global pass_fail_maxdiff;
    %Run the C demodulator
    cbits = fsk_demod_c(Fs,Rs,f1,f2,mod)';
    %Set up and run the octave demodulator
    states = fsk_horus_init(Fs,Rs);
    states.f1_tx = f1;
    states.f2_tx = f2;
    states.dA = 1;
    states.dF = 0;
    
    modin = mod;
    obits = [];
    while length(modin)>=states.nin
        ninold = states.nin;
        [bitbuf,states] = fsk_horus_demod(states, modin(1:states.nin));
        modin=modin(ninold+1:length(modin));
        obits = [obits bitbuf];
    end
    
    diff = sum(xor(obits,cbits))
    
    %Allow 2 bit difference between model and C
    pass = diff < 3
    if !pass
        printf('Demod failed test %s!\n',tname);
        figure(10);
        plot((1:length(obits)),obits,(1:length(cbits)),cbits);
        figure(11);
        plot(xor(obits,cbits));
    end
endfunction

% Random bit modulator test
% Pass random bits through the modulators and compare
function pass = test_mod_horuscfg_randbits
    bits = rand(1,10000)>.5;
    [dmod,cmod,omod,pass] = fsk_mod_test(8000,100,1200,1600,bits,"mod horuscfg randbits");
    
    figure(1)
    plot(dmod)
    title("Difference between octave and C mod impl");
    
    figure(2)
    %plot((1:length(cmod)),cmod,(1:length(omod)),omod,'mod random bits')
endfunction

% A big ol' channel impairment tester
% Shamlessly taken from fsk_horus
function pass = tfsk_run_sim(test_frame_mode,EbNodB,timing_offset,fading,df,dA)
  frames = 60;
  %EbNodB = 10;
  %timing_offset = 2.0; % see resample() for clock offset below
  %fading = 0;          % modulates tx power at 2Hz with 20dB fade depth, 
                       % to simulate balloon rotating at end of mission
  %df     = 0;          % tx tone freq drift in Hz/s
  %dA     = 1;          % amplitude imbalance of tones (note this affects Eb so not a gd idea)

  more off
  rand('state',1); 
  randn('state',1);

  % ----------------------------------------------------------------------

  % sm2000 config ------------------------
  %states = fsk_horus_init(96000, 1200);
  %states.f1_tx = 4000;
  %states.f2_tx = 5200;

  if test_frame_mode == 4
    % horus rtty config ---------------------
    states = fsk_horus_init(8000, 100);
    states.f1_tx = 1200;
    states.f2_tx = 1600;
    states.tx_bits_file = "horus_tx_bits_rtty.txt"; % Octave file of bits we FSK modulate
    
  end
                               
  if test_frame_mode == 5
    % horus binary config ---------------------
    states = fsk_horus_init(8000, 100);
    states.f1_tx = 1200;
    states.f2_tx = 1600;
    %%%states.tx_bits_file = "horus_tx_bits_binary.txt"; % Octave file of bits we FSK modulate
	states.tx_bits_file = "horus_payload_rtty.txt"
  end

  % ----------------------------------------------------------------------

  states.verbose = 0;%x1;
  N = states.N;
  P = states.P;
  Rs = states.Rs;
  nsym = states.nsym;
  Fs = states.Fs;
  states.df = df;
  states.dA = dA;

  EbNo = 10^(EbNodB/10);
  variance = states.Fs/(states.Rs*EbNo);

  % set up tx signal with payload bits based on test mode

  if test_frame_mode == 1
     % test frame of bits, which we repeat for convenience when BER testing
    test_frame = round(rand(1, states.nsym));
    tx_bits = [];
    for i=1:frames+1
      tx_bits = [tx_bits test_frame];
    end
  end
  if test_frame_mode == 2
    % random bits, just to make sure sync algs work on random data
    tx_bits = round(rand(1, states.nsym*(frames+1)));
  end
  if test_frame_mode == 3
    % ...10101... sequence
    tx_bits = zeros(1, states.nsym*(frames+1));
    tx_bits(1:2:length(tx_bits)) = 1;
  end
 
  if (test_frame_mode == 4) || (test_frame_mode == 5)

    % load up a horus msg from disk and modulate that

    test_frame = load(states.tx_bits_file);
    ltf = length(test_frame);
    ntest_frames = ceil((frames+1)*nsym/ltf);
    tx_bits = [];
    for i=1:ntest_frames
      tx_bits = [tx_bits test_frame];
    end
  end

  tx = fsk_horus_mod(states, tx_bits);

  tx = resample(tx, 1000, 1001); % simulated 1000ppm sample clock offset

  if fading
     ltx = length(tx);
     tx = tx .* (1.1 + cos(2*pi*2*(0:ltx-1)/Fs))'; % min amplitude 0.1, -20dB fade, max 3dB
  end

  noise = sqrt(variance)*randn(length(tx),1);
  rx    = tx + noise;
  %rx = real(rx);
  %b1 = fir2(100, [0 4000 5200 48000]/48000, [1 1 0.5 0.5]);
  %rx = filter(b1,1,rx);
  %[b a] = cheby2(6,40,[3000 6000]/(Fs/2));
  %rx = filter(b,a,rx);
  %rx = sign(rx);
  %rx(find (rx > 1)) = 1;
  %rx(find (rx < -1)) = -1;

  % dump simulated rx file
  ftx=fopen("fsk_horus_100bd_binary.raw","wb"); rxg = rx*1000; fwrite(ftx, rxg, "short"); fclose(ftx);

  timing_offset_samples = round(timing_offset*states.Ts);
  st = 1 + timing_offset_samples;
  rx_bits_buf = zeros(1,2*nsym);
  x_log = [];
  norm_rx_timing_log = [];
  f1_int_resample_log = [];
  f2_int_resample_log = [];
  f1_log = f2_log = [];
  EbNodB_log = [];
  rx_bits_log = [];
  rx_bits_sd_log = [];

  for f=1:frames

    % extract nin samples from input stream

    nin = states.nin;
    en = st + states.nin - 1;
    sf = rx(st:en);
    st += nin;

    % demodulate to stream of bits

    [rx_bits states] = fsk_horus_demod(states, sf);
    rx_bits_buf(1:nsym) = rx_bits_buf(nsym+1:2*nsym);
    rx_bits_buf(nsym+1:2*nsym) = rx_bits;
    rx_bits_log = [rx_bits_log rx_bits];
    rx_bits_sd_log = [rx_bits_sd_log states.rx_bits_sd];

    norm_rx_timing_log = [norm_rx_timing_log states.norm_rx_timing];
    x_log = [x_log states.x];
    f1_int_resample_log = [f1_int_resample_log abs(states.f1_int_resample)];
    f2_int_resample_log = [f2_int_resample_log abs(states.f2_int_resample)];
    f1_log = [f1_log states.f1];
    f2_log = [f2_log states.f2];
    EbNodB_log = [EbNodB_log states.EbNodB];

    if test_frame_mode == 1
       states = ber_counter(states, test_frame, rx_bits_buf);
    end
  end
  
  test_name = sprintf("tfsk_run_sim EbNodB:%d frames:%d timing_offset:%d fading:%d",EbNodB,frames,timing_offset,fading);
  [obits cbits pass]=fsk_demod_test(Fs,Rs,states.f1_tx,states.f2_tx,rx,test_name); 
  printf("Test %s done\n",test_name);
  
  if test_frame_mode == 1
    printf("frames: %d Tbits: %d Terrs: %d BER %4.3f\n", frames, states.Tbits,states. Terrs, states.Terrs/states.Tbits);
  end

  if test_frame_mode == 4
    extract_and_print_rtty_packets(states, rx_bits_log, rx_bits_sd_log)
  end

  if test_frame_mode == 5
    extract_and_decode_binary_packets(states, rx_bits_log);
  end

  figure(1);
  plot(f1_int_resample_log,'+')
  hold on;
  plot(f2_int_resample_log,'g+')
  hold off;

  figure(2)
  clf
  m = max(abs(x_log));
  plot(x_log,'+')
  axis([-m m -m m])
  title('fine timing metric')

  figure(3)
  clf
  subplot(211)
  plot(norm_rx_timing_log);
  axis([1 frames -1 1])
  title('norm fine timing')
  subplot(212)
  plot(states.nerr_log)
  title('num bit errors each frame')

  figure(4)
  clf
  subplot(211)
  plot(real(rx(1:Fs)))
  title('rx signal at demod input')
  subplot(212)
  plot(abs(fft(rx(1:Fs))))

  figure(5)
  clf
  plot(f1_log)
  hold on;
  plot(f2_log,'g');
  hold off;
  title('tone frequencies')
  axis([1 frames 0 Fs/2])

  figure(6)
  clf
  plot(EbNodB_log);
  title('Eb/No estimate')
  if(pass)
	close all;
  end
endfunction

function pass = ebno_battery_test(fading,df,dA)
    pass = 1
    ebnodbrange = fliplr(3:13)
    timing_offset=0;
  
    npass=zeros(1,length(ebnodbrange));
    for i=(1:length(npass))
        npass(i) = tfsk_run_sim(5,ebnodbrange(i),0,0,df,dA)
        assert(npass(i))
        pass = npass(i) && pass
    end
endfunction

function pass = test_fading_var(timing_offset,df,dA)
 
    pass = 1
    npass = ebno_battery_test(0,df,dA)
    assert(npass)
    pass = npass && pass
    npass = ebno_battery_test(1,df,dA)
    assert(npass)
    pass = npass && pass
endfunction


function pass = test_fsk_battery()
    pass = 1
    npass = test_mod_horuscfg_randbits
    pass = npass && pass
    assert(pass)
    npass = test_fading_var(0,0,1);
    pass = npass && pass
    assert(pass)
    if pass
        printf("***** All tests passed! *****\n");
    end
endfunction
