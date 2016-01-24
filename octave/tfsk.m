% tfsk.m
% Author: Brady O'Brien 8 January 2016



%   Copyright 2016 David Rowe
%  
%  All rights reserved.
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License version 2, as
%  published by the Free Software Foundation.  This program is
%  distributed in the hope that it will be useful, but WITHOUT ANY
%  WARRANTY; without even the implied warranty of MERCHANTABILITY or
%  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
%  License for more details.
%
%  You should have received a copy of the GNU Lesser General Public License
%  along with this program; if not, see <http://www.gnu.org/licenses/>.


% Octave script to check c port of fsk_horus against the fsk_horus.m
%
% [X] - Functions to wrap around fsk_mod and fsk_demod executables
%     [X] - fsk_mod
%     [X] - fsk_demod
% [X] - Functions to wrap around octave and c implementations, pass
%       same dataset, compare outputs, and give clear go/no-go
%     [X] - fsk_mod_test
%     [X] - fsk_demod_test
% [X] - Port of run_sim and EbNodB curve test battery
% [X] - Extract and compare more parameters from demod
% [X] - Run some tests in parallel

%
% FSK Modem test instructions --
% 1 - Compile tfsk.c and fm_mod.c
%     - tfsk.c is in unittest/, so build must not be configured for release
% 2 - Change tfsk_location and fsk_mod_location to point to tfsk
% 3 - Ensure octave packages signal and parallel are installed
% 4 - run tfsk.m. It will take care of the rest.
%

%tfsk executable path/file
global tfsk_location = '../build/unittest/tfsk';



fsk_horus_as_a_lib = 1; % make sure calls to test functions at bottom are disabled
fsk_horus_2fsk;  
pkg load signal;
pkg load parallel;
graphics_toolkit('gnuplot');


global mod_pass_fail_maxdiff = 1e-3/50000;

function mod = fsk_mod_c(Fs,Rs,f1,f2,bits)
    global tfsk_location;
    %command to be run by system to launch the modulator
    command = sprintf('%s M %d %d %d %d fsk_mod_ut_bitvec fsk_mod_ut_modvec fsk_mod_ut_log.txt',tfsk_location,f1,f2,Fs,Rs);
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


%Compare 2 vectors, fail if they are not close enough
function pass = vcompare(vc,voct,vname,tname,tol)
    
    %Get delta of vectors
    dvec = abs(abs(vc)-abs(voct));     
    
    %Normalize difference
    dvec = dvec ./ abs(max(abs(voct)));
    
    maxdvec = abs(max(dvec));
    pass = maxdvec<tol;
    
    printf('Comparing vectors %s in test %s. Diff is %f\n',vname,tname,maxdvec);
    
    if pass == 0
        printf('\n*** vcompare failed %s in test %s. Diff: %f Tol: %f\n\n',vname,tname,maxdvec,tol);
        
        titlestr = sprintf('Diff between C and Octave of %s for %s',vname,tname)
        figure(12)
        plot(abs(dvec))
        title(titlestr)
        
        figure(13)
        plot(abs(va))
    
        figure(14)
        plot(abs(vb))
    end
    assert(pass);
    
endfunction

function test_stats = fsk_demod_xt(Fs,Rs,f1,f2,mod,tname)
    global tfsk_location;
    %Name of executable containing the modulator
    fsk_demod_ex_file = '../build/unittest/tfsk';
    modvecfilename = sprintf('fsk_demod_ut_modvec_%d',getpid());
    bitvecfilename = sprintf('fsk_demod_ut_bitvec_%d',getpid());
    tvecfilename = sprintf('fsk_demod_ut_tracevec_%d.txt',getpid());
    
    %command to be run by system to launch the demod
    command = sprintf('%s D %d %d %d %d %s %s %s',tfsk_location,f1,f2,Fs,Rs,modvecfilename,bitvecfilename,tvecfilename);
    
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
    o_ppm = [];
    o_rx_timing = [];
    %Run octave demod, dump some test vectors
    states = fsk_horus_init(Fs,Rs);
    Ts = states.Ts;
    P = states.P;
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
        o_f1_dc = [o_f1_dc states.f1_dc(1:states.Nmem-Ts/P)];
        o_f2_dc = [o_f2_dc states.f2_dc(1:states.Nmem-Ts/P)];
        o_f1_int = [o_f1_int states.f1_int];
        o_f2_int = [o_f2_int states.f2_int];
        o_EbNodB = [o_EbNodB states.EbNodB];
        o_ppm = [o_ppm states.ppm];
        o_rx_timing = [o_rx_timing states.rx_timing];
    end
    
    % One part-per-thousand allowed on important parameters
    pass =         vcompare(o_f1_dc,      t_f1_dc,    'f1_dc',    tname,.001);
    pass = pass && vcompare(o_f2_dc,      t_f2_dc,    'f2_dc',    tname,.001);
    pass = pass && vcompare(o_f1_int,     t_f1_int,   'f1_int',   tname,.001);
    pass = pass && vcompare(o_f2_int,     t_f2_int,   'f2_int',   tname,.001);
    pass = pass && vcompare(o_rx_timing,  t_rx_timing,'rx_timing',tname,.001);
    
    % Much larger tolerances on unimportant statistics
    pass = pass && vcompare(o_EbNodB,     t_EbNodB,   'EbNodB',   tname,.05);
    pass = pass && vcompare(o_ppm   ,     t_ppm,      'ppm',      tname,.02);
    
    diffpass = sum(xor(obits,bits'))<4;
    diffbits = sum(xor(obits,bits'));
    
    if diffpass==0
        printf('\n***bitcompare test failed test %s diff %d\n\n',tname,sum(xor(obits,bits')))
        figure(15)
        plot(xor(obits,bits'))
        title(sprintf('Bitcompare failure test %s',tname))
    end
    
    pass = pass && diffpass;
    
    
    test_stats.pass = pass;
    test_stats.diff = sum(xor(obits,bits'));
    test_stats.cbits = bits';
    test_stats.obits = obits;
    
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
    omod = fsk_horus_mod(states,bits');
    
    dmod = cmod-omod;
    pass = max(dmod)<(mod_pass_fail_maxdiff*length(dmod))
    if !pass
        printf('Mod failed test %s!\n',tname);
    end
endfunction

% Random bit modulator test
% Pass random bits through the modulators and compare
function pass = test_mod_horuscfg_randbits
    rand('state',1); 
    randn('state',1);
    bits = rand(1,10000)>.5;
    [dmod,cmod,omod,pass] = fsk_mod_test(8000,100,1200,1600,bits,"mod horuscfg randbits");
    
    if(!pass)
        figure(1)
        plot(dmod)
        title("Difference between octave and C mod impl");
    end
    
endfunction

% A big ol' channel impairment tester
% Shamlessly taken from fsk_horus
% This throws some channel imparment or another at the C and octave modem so they 
% may be compared.
function stats = tfsk_run_sim(test_frame_mode,EbNodB,timing_offset,fading,df,dA)
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
	states.tx_bits_file = "horus_payload_rtty.txt";
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
  
  test_name = sprintf("tfsk_run_sim EbNodB:%d frames:%d timing_offset:%d fading:%d df:%d",EbNodB,frames,timing_offset,fading,df);
  tstats = fsk_demod_xt(Fs,Rs,states.f1_tx,states.f2_tx,rx,test_name); 
  printf("Test %s done\n",test_name);
  
  pass = tstats.pass
  obits = tstats.obits;
  cbits = tstats.cbits;
  
  % Figure out BER of octave and C modems
  bitcnt = length(tx_bits)
  rx_bits = obits;
  ber = 1;
  ox = 1;
  for offset = (1:100)
    nerr = sum(xor(rx_bits(offset:length(rx_bits)),tx_bits(1:length(rx_bits)+1-offset)));
    bern = nerr/(bitcnt-offset);
    if(bern < ber)
      ox = offset;
      best_nerr = nerr;
    end
    ber = min([ber bern]);
  end
  offset = ox;
  bero = ber;
  ber = 1;
  rx_bits = cbits;
  ox = 1;
  for offset = (1:100)
    nerr = sum(xor(rx_bits(offset:length(rx_bits)),tx_bits(1:length(rx_bits)+1-offset)));
    bern = nerr/(bitcnt-offset);
    if(bern < ber)
      ox = offset;
      best_nerr = nerr;
    end
    ber = min([ber bern]);
  end
  offset = ox;
  berc = ber;
  stats.berc = berc;
  stats.bero = bero;
  
  % non-coherent BER theory calculation
  % It was complicated, so I broke it up

  ms = 2;
  ns = (1:ms-1);
  as = (-1).^(ns+1);
  bs = (as./(ns+1));
  
  cs = ((ms-1)./ns);

  ds = ns.*log2(ms);
  es = ns+1;
  fs = exp( -(ds./es)*EbNo );
  
  thrncoh = ((ms/2)/(ms-1)) * sum(bs.*((ms-1)./ns).*exp( -(ds./es)*EbNo ));
  
  stats.thrncoh = thrncoh;
  stats.pass = pass;
endfunction


function pass = ebno_battery_test(timing_offset,fading,df,dA)
    %Range of EbNodB over which to test
    ebnodbrange = (5:13);
    ebnodbs = length(ebnodbrange);
    
    mode = 5;
    %Replication of other parameters for parcellfun
    modev   = repmat(mode,1,ebnodbs);
    timingv = repmat(timing_offset,1,ebnodbs);
    fadingv = repmat(fading,1,ebnodbs);
    dfv     = repmat(df,1,ebnodbs);
    dav     = repmat(dA,1,ebnodbs);
    
    statv = pararrayfun(floor(1.25*nproc()),@tfsk_run_sim,modev,ebnodbrange,timingv,fadingv,dfv,dav);
    %statv = arrayfun(@tfsk_run_sim,modev,ebnodbrange,timingv,fadingv,dfv,dav);

    passv = zeros(1,length(statv));
    for ii=(1:length(statv))
        passv(ii)=statv(ii).pass;
    end
    
    %All pass flags are '1'
    pass = sum(passv)>=length(passv);
    %and no tests died
    pass = pass && length(passv)==ebnodbs;
    passv
    assert(pass)
endfunction

%Test with and without channel fading
function pass = test_fading_var(timing_offset,df,dA)
    pass = ebno_battery_test(timing_offset,1,df,dA)
    assert(pass)
    pass = pass && ebno_battery_test(timing_offset,0,df,dA)
    assert(pass)
endfunction

%Test with and without sample clock offset
function pass = test_timing_var(df,dA)
    pass = test_fading_var(1,df,dA)
    assert(pass)
    pass = pass && test_fading_var(0,df,dA)
    assert(pass)
endfunction

%Test with and without 1 Hz/S freq drift
function pass = test_drift_var()
    pass = test_timing_var(1,1)
    assert(pass)
    pass = pass && test_timing_var(0,1)
    assert(pass)
endfunction

function pass = test_fsk_battery()
    pass = test_mod_horuscfg_randbits
    assert(pass)
    pass = pass && test_drift_var();
    assert(pass)
    
    if pass
        printf("***** All tests passed! *****\n");
    end
endfunction

function plot_fsk_bers()
    %Range of EbNodB over which to plot
    ebnodbrange = (5:13);
    
    berc = ones(1,length(ebnodbrange));
    bero = ones(1,length(ebnodbrange));
    berinc = ones(1,length(ebnodbrange));
    ebnodbs = length(ebnodbrange)
    mode = 5;
    %Replication of other parameters for parcellfun
    modev   = repmat(mode,1,ebnodbs);
    timingv = repmat(0,1,ebnodbs);
    fadingv = repmat(0,1,ebnodbs);
    dfv     = repmat(0,1,ebnodbs);
    dav     = repmat(1,1,ebnodbs);
    statv = pararrayfun(floor(1.25*nproc()),@tfsk_run_sim,modev,ebnodbrange,timingv,fadingv,dfv,dav);
    
    for ii = (1:length(statv))
        stat = statv(ii);
        berc(ii)=stat.berc;
        bero(ii)=stat.bero;
        berinc(ii)=stat.thrncoh;
    end
    
    close all
    figure(2);
    clf;
    semilogy(ebnodbrange, berinc,'r;2FSK non-coherent theory;')
    hold on;
    semilogy(ebnodbrange, bero ,'g;Octave fsk horus 2FSK Demod;')
    semilogy(ebnodbrange, berc, 'v;C fsk horus 2FSK Demod;')
    hold off;
    grid("minor");
    axis([min(ebnodbrange) max(ebnodbrange) 1E-5 1])
    legend("boxoff");
    xlabel("Eb/No (dB)");
    ylabel("Bit Error Rate (BER)")
 
endfunction

plot_fsk_bers
test_fsk_battery
