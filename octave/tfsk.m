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
pkg load signal;

graphics_toolkit('gnuplot');


global mod_pass_fail_maxdiff = 1e-3/50000;
global demod_pass_fail_maxdiff = .1;

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
    fsk_demod_ex_file = '../build/unittest/tfsk';
    modvecfilename = sprintf('fsk_demod_ut_modvec_%d',getpid());
    bitvecfilename = sprintf('fsk_demod_ut_bitvec_%d',getpid());
    tvecfilename = sprintf('fsk_demod_ut_tracevec_%d.txt',getpid());
    
    %command to be run by system to launch the demod
    command = sprintf('%s %d %d %d %d %s %s %s',fsk_demod_ex_file,Fs,Rs,f1,f2,modvecfilename,bitvecfilename,tvecfilename);
    
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
    %delete(tvecfilename);
endfunction

%Compare 2 vectors, fail if they are not close enough
function pass = vcompare(va,vb,vname,tname)
    global demod_pass_fail_maxdiff;
    
    %Get delta of vectors
    dvec = abs(va)-abs(vb);     
    
    %Normalize difference
    dvec = dvec ./ max(va);
    
    titlestr = sprintf('Diff between C and Octave of %s for %s',vname,tname);
    pass = max(dvec)<(demod_pass_fail_maxdiff)
    maxdvec = max(dvec)
    %figure(12)
    %title(titlestr)
    %plot(abs(dvec))
    
    %figure(13)
    %plot(abs(va))
    
    %figure(14)
    %plot(abs(vb))
    
    printf('Comparing vectors %s in test %s\n',vname,tname);
    assert(pass);
    
endfunction

function test_stats = fsk_demod_xt(Fs,Rs,f1,f2,mod,tname)
    global mod_pass_fail_maxdiff;
    %Name of executable containing the modulator
    fsk_demod_ex_file = '../build/unittest/tfsk';
    modvecfilename = sprintf('fsk_demod_ut_modvec_%d',getpid());
    bitvecfilename = sprintf('fsk_demod_ut_bitvec_%d',getpid());
    tvecfilename = sprintf('fsk_demod_ut_tracevec_%d.txt',getpid());
    
    %command to be run by system to launch the demod
    command = sprintf('%s %d %d %d %d %s %s %s',fsk_demod_ex_file,Fs,Rs,f1,f2,modvecfilename,bitvecfilename,tvecfilename);
    
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
    
    %Load test vec dump
    load(tvecfilename);
    
    %Clean up files
    delete(bitvecfilename);
    delete(modvecfilename);
    delete(tvecfilename);
    
    
    o_f1_dc = [];
    o_f2_dc = [];
    o_f1_int = [];
    o_f2_int = [];
    o_EbNodB = [];
    o_rx_timing = [];
    %Run octave demod, dump some test vectors
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
        
        %Save other parameters
        o_f1_dc = [o_f1_dc states.f1_dc];
        o_f2_dc = [o_f2_dc states.f2_dc];
        o_f1_int = [o_f1_int states.f1_int];
        o_f2_int = [o_f2_int states.f2_int];
        o_EbNodB = [o_EbNodB states.EbNodB];
        o_rx_timing = [o_rx_timing states.rx_timing];
    end
    
    pass = 1;
    pass = and(pass,vcompare(o_rx_timing,t_rx_timing,'rx_timing',tname));
    pass = and(pass,vcompare(o_f1_int,t_f1_int,'f1_int',tname));
    pass = and(pass,vcompare(o_f2_int,t_f2_int,'f2_int',tname));
     
    diff = sum(xor(obits,bits'))
    
    pass = pass && diff<2;
    test_stats.pass = pass;
    assert(pass)
    test_stats.diff = diff;
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
    rand('state',1); 
    randn('state',1);
    bits = rand(1,10000)>.5;
    [dmod,cmod,omod,pass] = fsk_mod_test(8000,100,1200,1600,bits,"mod horuscfg randbits");
    
    figure(1)
    plot(dmod)
    title("Difference between octave and C mod impl");
    
endfunction

% A big ol' channel impairment tester
% Shamlessly taken from fsk_horus
% This throws some channel imparment or another at the C and octave modem so they 
% may be compared.
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
  
  if test_frame_mode == 2
    % horus rtty config ---------------------
    states = fsk_horus_init(8000, 100);
    states.f1_tx = 1200;
    states.f2_tx = 1600;
    
  end

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

  if timing_offset
    tx = resample(tx, 1000, 1001); % simulated 1000ppm sample clock offset
  end
  
  if fading
     ltx = length(tx);
     tx = tx .* (1.1 + cos(2*pi*2*(0:ltx-1)/Fs))'; % min amplitude 0.1, -20dB fade, max 3dB
  end

  noise = sqrt(variance)*randn(length(tx),1);
  rx    = tx + noise;
  
  test_name = sprintf("tfsk_run_sim EbNodB:%d frames:%d timing_offset:%d fading:%d",EbNodB,frames,timing_offset,fading);
  tstats = fsk_demod_xt(Fs,Rs,states.f1_tx,states.f2_tx,rx,test_name); 
  printf("Test %s done\n",test_name);
  
  pass = tstats.pass
  
endfunction

function pass = ebno_battery_test(fading,df,dA)
    pass = 1
    ebnodbrange = (2:13)
    timing_offset=0;
  
    npass=zeros(1,length(ebnodbrange));
    for i=(1:length(npass))
        npass(i) = tfsk_run_sim(5,ebnodbrange(i),0,0,df,dA)
        assert(npass(i))
        pass = npass(i) && pass
    end
endfunction

function pass = ebno_battery_test(timing_offset,fading,df,dA)
    %Range of EbNodB over which to test
    ebnodbrange = fliplr(3:13);
    ebnodbs = length(ebnodbrange);
    
    mode = 2;
    %Replication of other parameters for parcellfun
    modev   = repmat(mode,1,ebnodbs);
    timingv = repmat(timing_offset,1,ebnodbs);
    fadingv = repmat(fading,1,ebnodbs);
    dfv     = repmat(df,1,ebnodbs);
    dav     = repmat(dA,1,ebnodbs);
    
    passv = pararrayfun(floor(.7*nproc()),@tfsk_run_sim,modev,ebnodbrange,timingv,fadingv,dfv,dav);
    %passv = arrayfun(@tfsk_run_sim,modev,ebnodbrange,timingv,fadingv,dfv,dav);
    
    pass = sum(passv)>=length(passv)
    passv
    assert(pass)
endfunction

%Test with and without channel fading
function pass = test_fading_var(timing_offset,df,dA)
    pass = ebno_battery_test(timing_offset,0,df,dA)
    assert(pass)
    pass = pass && ebno_battery_test(timing_offset,1,df,dA)
    assert(pass)
endfunction

%Test with and without sample clock offset
function pass = test_timing_var(df,dA)
    pass = test_fading_var(0,df,dA)
    assert(pass)
    pass = pass && test_fading_var(1,df,dA)
    assert(pass)
endfunction

function pass = test_fsk_battery()
    pass = test_mod_horuscfg_randbits
    assert(pass)
    pass = pass && test_timing_var(0,1);
    assert(pass)
    
    if pass
        close all;
        printf("***** All tests passed! *****\n");
    end
endfunction
