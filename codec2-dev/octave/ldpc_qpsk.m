% ldpc_qpsk.m
%
% David Rowe 18 Dec 2013
%
% Similar to ldpc_short.m, but derived from ldpcut.m and uses QPSK and
% CML 2D functunctions and QPSK.  Probably should combine this and
% ldpc_short.m some day.

% Our LDPC library

ldpc;
qpsk;
gp_interleaver;


function sim_out = run_simulation(sim_in)

  % Note this is effective Eb/No of payload data bits, sorta thing we
  % plot on BER versus Eb/No graphs of decoded data.  So if we have a
  % rate 1/2 code, each codeword bit will have Eb/No - 3dB.

  EbNodBvec = sim_in.EbNodBvec;

  % for wimax code frame size specifies code

  if isfield(sim_in, "framesize")
    framesize = sim_in.framesize;
    rate = sim_in.rate; 
  end

  Ntrials = sim_in.Ntrials;
  verbose = sim_in.verbose;
  if isfield(sim_in, "hf_en")
    hf_en = sim_in.hf_en;
  else
    hf_en = 0;
  end
  ldpc_code = sim_in.ldpc_code;
  interleave_en = sim_in.interleave_en;

  % Init LDPC code ------------------------------------

  mod_order = 4; bps = 2;
  modulation = 'QPSK';
  mapping = 'gray';

  demod_type = 0;
  decoder_type = 0;
  max_iterations = 100;

  if ldpc_code == 1
    code_param = ldpc_init_wimax(rate, framesize, modulation, mod_order, mapping);
  end
  if ldpc_code == 0
    load HRA_112_112.txt
    [code_param framesize rate] = ldpc_init_user(HRA_112_112, modulation, mod_order, mapping);
  end
  if ldpc_code == 2
    load('H2064_516_sparse.mat');
    HRA = full(HRA);  
    [code_param framesize rate] = ldpc_init_user(HRA, modulation, mod_order, mapping);
  end 
  if ldpc_code == 3
    load('h0p25d.mat');
    %HRA = full(HRA);  
    [code_param framesize rate] = ldpc_init_user(H, modulation, mod_order, mapping);
  end

  % set up optional HF (multipath) model ------------------------------------

  % signal is arranged as Nc parallel carriers.  Nc is chosen such
  % that payload data rate is 700 bits/s.  So for higher rate codes Nc
  % will be smaller.

  Rs = 50;
  vocoder_bps = 700; raw_bps = vocoder_bps/rate;  
  Nc = round(raw_bps/(Rs*bps));
  Tp = (framesize/Nc)/Rs; Tp_codec2 = 0.04;
  fading = ones(1,Ntrials*code_param.code_bits_per_frame/bps);

  printf("framesize: %d  rate: %3.2f  Nc: %d\n", framesize, rate, Nc);

  if hf_en

    % We assume symbols spread acroos Nc OFDM carriers

    dopplerSpreadHz = 1.0; path_delay = 1E-3*Rs;

    if isfield(sim_in, "dopplerSpreadHz") 
      dopplerSpreadHz = sim_in.dopplerSpreadHz;
    end
    if isfield(sim_in, "path_delay") 
      path_delay = sim_in.path_delay;
    end
    printf("Doppler Spread: %3.2f Hz Path Delay: %3.2f symbols\n", dopplerSpreadHz, path_delay);

    % reset seed so we get same fading channel on every simulation

    randn('seed',1);

    Ns = Ntrials*code_param.code_bits_per_frame/bps;
    Nr = ceil(Ns/Nc);
    hf_model = zeros(Nr,Nc); 

    % note we ask for 10% more samples than we need, as
    % doppler_spread() function can sometimes return slightly less
    % than we need due to round off

    spread1 = doppler_spread(dopplerSpreadHz, Rs, Nr*1.1);
    spread2 = doppler_spread(dopplerSpreadHz, Rs, Nr*1.1);
    spread1 = spread1(1:Nr); 
    spread2 = spread2(1:Nr);
    hf_gain = 1.0/sqrt(var(spread1)+var(spread2));
  end

  % ----------------------------------
  % run simulation at each Eb/No point
  % ----------------------------------

  for ne = 1:length(EbNodBvec)
    randn('seed',1);
    rand('seed',1);

    % Given Eb/No of payload data bits, work out Es/No we need to
    % apply to each channel symbol:
    %
    % i) Each codeword bit gets noise: Eb/No - 3 (for a rate 1/2 code)
    % ii) QPSK means two bits/symbol.: Es/No = Eb/No + 3
    %
    % -> which neatly cancel out ...... (at least for rate 1/2)

    EsNodB = EbNodBvec(ne) + 10*log10(rate) + 10*log10(bps);
    EsNo = 10^(EsNodB/10);
    variance = 1/EsNo;
    hf_r = 1;
    
    Tbits = Terrs = Ferrs = Terrs_raw = Tbits_raw = 0;
    
    tx_bits = [];
    tx_symbols = []; 
    rx_symbols = [];

    % Encode a bunch of frames

    for nn=1:Ntrials        
      atx_bits = round(rand( 1, code_param.data_bits_per_frame));
      tx_bits = [tx_bits atx_bits];
      [tx_codeword atx_symbols] = ldpc_enc(atx_bits, code_param);
      if interleave_en
        atx_symbols = gp_interleave(atx_symbols);
      end
      tx_symbols = [tx_symbols atx_symbols];
    end
    
    rx_symbols = tx_symbols;

    % Optional HF (multipath) channel model

    if hf_en

      % Simplified rate Rs OFDM simulation model that doesn't
      % include ISI, just freq filtering.  We assume perfect phase
      % estimation so it's just amplitude distortion.  We distribute
      % symbols across Nc carriers

      Ns = length(tx_symbols);
      w = (1:Nc)*2*pi;  
      rx_symbols = [rx_symbols zeros(1,Nr*Nc-Ns)];  % pad out to integer number of rows

      for r=1:Nr
        for c=1:Nc
          hf_model(hf_r,c) = hf_gain*(spread1(hf_r) + exp(-j*w(c)*path_delay)*spread2(hf_r));
          rx_symbols(Nc*(r-1)+c) *= abs(hf_model(hf_r,c));
          fading(Nc*(r-1)+c) = abs(hf_model(hf_r,c));
        end
        hf_r++;
      end
      rx_symbols = rx_symbols(1:Ns); % remove padding
    end

    % Add AWGN noise, 0.5 factor splits power evenly between Re & Im

    noise = sqrt(variance*0.5)*(randn(1,length(tx_symbols)) + j*randn(1,length(tx_symbols)));
    rx_symbols += noise;

    % Decode a bunch of frames

    rx_bits_log = []; error_positions_log = [];
    Nerrs_raw_log = [];

    for nn = 1: Ntrials        
      st = (nn-1)*code_param.symbols_per_frame + 1;
      en = (nn)*code_param.symbols_per_frame;

      % coded 

      arx_symbols = rx_symbols(st:en);
      afading = fading(st:en);
      if interleave_en
        arx_symbols = gp_deinterleave(arx_symbols);
        afading = gp_deinterleave(afading);
      end

      arx_codeword = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, arx_symbols, EsNo, afading);
      st = (nn-1)*code_param.data_bits_per_frame + 1;
      en = (nn)*code_param.data_bits_per_frame;
      arx_bits = arx_codeword(1:code_param.data_bits_per_frame);
      error_positions = xor(arx_bits, tx_bits(st:en));
      error_positions_log = [error_positions_log error_positions];
      Nerrs = sum( error_positions);
      rx_bits_log = [rx_bits_log arx_bits];
        
      % uncoded - to est raw BER compare first half or received frame to tx_bits as code is systematic
      
      raw_rx_bits = [];
      for s=1:code_param.symbols_per_frame*rate
        raw_rx_bits = [raw_rx_bits qpsk_demod(arx_symbols(s))];
      end
      Nerrs_raw = sum(xor(raw_rx_bits, tx_bits(st:en)));
      Nbits_raw = code_param.data_bits_per_frame;
      Nerrs_raw_log = [Nerrs_raw_log Nerrs_raw];

      if verbose == 2
        % print "." if frame decoded without errors, 'x' if we can't decode

        if Nerrs > 0, printf('x'),  else printf('.'),  end
        if (rem(nn, 50)==0),  printf('\n'),  end    
      end

      if Nerrs > 0,  Ferrs = Ferrs + 1;  end
      Terrs     += Nerrs;
      Tbits     += code_param.data_bits_per_frame;        
      Terrs_raw += Nerrs_raw;
      Tbits_raw += Nbits_raw;
    end

    % Alternative Codec 2 packet rate measurement indep of framesize

    Nerrs_codec2_log = []; Ncodecpacketsize = 28;
    Perrs = 0; Npackets = floor(length(tx_bits)/Ncodecpacketsize);
    for p=1:Ncodecpacketsize:Npackets*Ncodecpacketsize
      Nerrs = sum(xor(tx_bits(p:p+Ncodecpacketsize-1), rx_bits_log(p:p+Ncodecpacketsize-1)));
      if Nerrs
        Perrs++;
      end
      Nerrs_codec2_log = [Nerrs_codec2_log Nerrs];
    end

    if verbose
      printf("\nCoded EbNodB: %3.2f BER: %4.3f Tbits: %6d Terrs: %6d FER: %4.3f Tframes: %d Ferrs: %d\n",
             EbNodBvec(ne), Terrs/Tbits, Tbits, Terrs,  Ferrs/Ntrials, Ntrials, Ferrs);
      EbNodB_raw = EbNodBvec(ne) + 10*log10(rate);
      printf("Raw EbNodB..: %3.2f BER: %4.3f Tbits: %6d Terrs: %6d\n", 
             EbNodB_raw, Terrs_raw/Tbits_raw, Tbits_raw, Terrs_raw);
      printf("Codec 2 PER: %5.4f Npackets: %d Perrs: %d\n", Perrs/Npackets, Npackets, Perrs);
    end

    sim_out.rate = rate;
    sim_out.BER(ne) = Terrs/Tbits;
    sim_out.PER(ne) = Perrs/Npackets;
    sim_out.error_positions = error_positions_log;

    if hf_en && (verbose > 1)
      figure(2); clf;
      plot(real(rx_symbols(Ns/2:Ns)), imag(rx_symbols(Ns/2:Ns)), '+');
      axis([-2 2 -2 2]);
      title('Scatter')

      figure(3); clf;
      subplot(211);
      stem((1:Ntrials)*Tp, Nerrs_raw_log);
      subplot(212);
      stem((1:Npackets)*Tp_codec2, Nerrs_codec2_log);

      figure(4); clf;

      % limit mesh plot to Np points to plot quickly
      
      Np = 500;
      step = ceil(hf_r/Np);
      mesh(1:Nc, (1:step:hf_r-1)/Rs, abs(hf_model(1:step:hf_r-1,:)))
      title('HF channel amplitude');
      xlabel('Carrier');
      ylabel('Time (s)');

      figure(5)
      subplot(211); plot((1:hf_r-1)/Rs, abs(spread1(1:hf_r-1)));
      subplot(212); plot((1:hf_r-1)/Rs, abs(spread2(1:hf_r-1)));
      title('HF spreading function amplitudes')
    end
  end

endfunction


% ---------------------------------------------------------------------------------
% Run a bunch of trials at just one EbNo point
% ---------------------------------------------------------------------------------

function run_single(Nbits=700*10, EbNodB=9, hf_en=0, ldpc_code=1, framesize=576, interleave_en=0, error_pattern_filename)
  sim_in.ldpc_code = ldpc_code;

  if sim_in.ldpc_code == 0
    % Our HRA short LDPC code
    sim_in.rate=0.5;
    sim_in.framesize=448*4+448; 
  end
  if sim_in.ldpc_code == 1
    % CML wimax codes
    sim_in.rate = 0.5; 
    sim_in.framesize = framesize;
  end
  if sim_in.ldpc_code == 2
    sim_in.rate=0.8;
    sim_in.framesize=2064+516; 
  end 
  if sim_in.ldpc_code == 3
    sim_in.rate=0.25;
    sim_in.framesize=2064+516; 
  end 

  sim_in.verbose = 2;
  sim_in.Ntrials = ceil(Nbits/(sim_in.framesize*sim_in.rate));
  sim_in.EbNodBvec = EbNodB;
  sim_in.hf_en = hf_en;
  sim_in.interleave_en = interleave_en;

  sim_out = run_simulation(sim_in);

  if nargin == 7
    fep = fopen(error_pattern_filename, "wb");
    fwrite(fep, sim_out.error_positions, "short");
    fclose(fep);
  end

end

% ---------------------------------------------------------------------------------
% Lets draw some Eb/No versus BER curves 
% ---------------------------------------------------------------------------------

function plot_curves(Nbits=700*60)
  sim_in.EbNodBvec = -2:12;
  sim_in.verbose = 2;
  sim_in.interleave_en = 1;

  % Low rate 0.25 VK5DSP code

  sim_in.ldpc_code = 3;
  sim_in.rate = 0.25;
  sim_in.framesize = 448*4+448;  
  sim_in.Ntrials = floor(Nbits/(sim_in.framesize*sim_in.rate));

  sim_in.hf_en = 0; sim_out_awgn_low = run_simulation(sim_in);
  sim_in.hf_en = 1; sim_out_hf_low = run_simulation(sim_in);

  % Wimax code

  sim_in.ldpc_code = 1;
  sim_in.rate = 0.5; 
  sim_in.framesize = 576*4;
  sim_in.Ntrials = floor(Nbits/(sim_in.framesize*sim_in.rate));

  sim_in.hf_en = 0; sim_out_awgn_wimax = run_simulation(sim_in);
  sim_in.hf_en = 1; sim_out_hf_wimax = run_simulation(sim_in);

  % Our short code from VK5DSP

  sim_in.ldpc_code = 0;
  sim_in.rate = 0.5;
  sim_in.framesize = 224;  
  sim_in.Ntrials = floor(Nbits/(sim_in.framesize*sim_in.rate));

  sim_in.hf_en = 0; sim_out_awgn_short = run_simulation(sim_in);
  sim_in.hf_en = 1; sim_out_hf_short = run_simulation(sim_in);

  % Rate 0.8 Wenet code from VK5DSP

  sim_in.ldpc_code = 2;
  sim_in.rate = 0.8;
  sim_in.framesize = 2064+512;  
  sim_in.Ntrials = floor(Nbits/(sim_in.framesize*sim_in.rate));

  sim_in.hf_en = 0; sim_out_awgn_wenet = run_simulation(sim_in);
  sim_in.hf_en = 1; sim_out_hf_wenet = run_simulation(sim_in);

  % plots -------------------------

  EbNodB = sim_in.EbNodBvec;
  uncoded_awgn_ber_theory = 0.5*erfc(sqrt(10.^(EbNodB/10)));

  EbNoLin = 10.^(EbNodB/10);
  uncoded_hf_ber_theory = 0.5.*(1-sqrt(EbNoLin./(EbNoLin+1)));

  figure(1)
  clf
  semilogy(EbNodB, uncoded_awgn_ber_theory,'r-+;AWGN Uncoded;','markersize', 10, 'linewidth', 2)
  hold on;
  semilogy(EbNodB, uncoded_hf_ber_theory,'r-o;HF Uncoded;','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_awgn_wimax.BER+1E-10,'g-+;AWGN LDPC (2304,1152);','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_hf_wimax.BER+1E-10,'g-o;HF LDPC (2304,1152);','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_awgn_short.BER+1E-10,'b-+;AWGN LDPC (224,112);','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_hf_short.BER+1E-10,'b-o;HF LDPC (224,112);','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_awgn_wenet.BER+1E-10,'c-+;AWGN LDPC (2576,2064);','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_hf_wenet.BER+1E-10,'c-o;HF LDPC (2576,2064);','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_awgn_low.BER+1E-10,'k-+;AWGN LDPC (1792,448);','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_hf_low.BER+1E-10,'k-o;HF LDPC (1792,448);','markersize', 10, 'linewidth', 2);
  hold off;
  grid('minor')
  xlabel('Eb/No (dB)')
  ylabel('BER')
  axis([min(EbNodB) max(EbNodB) 1E-3 5E-1])
  legend('boxoff')
  epsname = sprintf("ldpc_qpsk_ber.eps");
  print('-deps', '-color', epsname)

  uncoded_awgn_per_theory = 1 - (1-uncoded_awgn_ber_theory).^28;
  uncoded_hf_per_theory = 1 - (1-uncoded_hf_ber_theory).^28;

  figure(2)
  clf
  semilogy(EbNodB, uncoded_awgn_per_theory,'r-+;AWGN Uncoded;','markersize', 10, 'linewidth', 2)
  hold on;
  semilogy(EbNodB, uncoded_hf_per_theory,'r-o;HF Uncoded;','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_awgn_wimax.PER+1E-10,'g-+;AWGN LDPC (2304,1152);','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_hf_wimax.PER+1E-10,'g-o;HF LDPC (2304,1152);','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_awgn_short.PER+1E-10,'b-+;AWGN LDPC (224,112);','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_hf_short.PER+1E-10,'b-o;HF LDPC (224,112);','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_awgn_wenet.PER+1E-10,'c-+;AWGN LDPC (2576,2064);','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_hf_wenet.PER+1E-10,'c-o;HF LDPC (2576,2064);','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_awgn_low.PER+1E-10,'k-+;AWGN LDPC (1792,448);','markersize', 10, 'linewidth', 2);
  semilogy(EbNodB, sim_out_hf_low.PER+1E-10,'k-o;HF LDPC (1792,448);','markersize', 10, 'linewidth', 2);
   hold off;
  grid('minor')
  xlabel('Eb/No (dB)')
  ylabel('PER')
  axis([min(EbNodB) max(EbNodB) 1E-2 1])
  legend('boxoff')
  legend("location", "southwest");
  epsname = sprintf("ldpc_qpsk_per.eps");
  print('-deps', '-color', epsname)
end


% ---------------------------------------------------------------------------------
% Start simulations here
% ---------------------------------------------------------------------------------

more off;
format;

init_cml('~/cml/');

%run_single(Nbits=700*5, EbNo=6, hf_en=1, ldpc_code=3, framesize=576*4, 1)
plot_curves(700*60);



