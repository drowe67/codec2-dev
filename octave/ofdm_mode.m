% ofdm_mode.m
%
% Library of functions to help setting up OFDM modes

%------------------------------------------------------------------------------
% ofdm_init_mode - Helper function to set up modems for various FreeDV modes,
%                  and parse mode string.
%------------------------------------------------------------------------------

1;

function config = ofdm_init_mode(mode="700D")
  % defaults for 700D
  
  Tcp = 0.002; 
  Ns = 8;
  Ts = 0.018; 
  Nc = 17;
  config.bps = 2; 
  config.Np = 1;
  config.Ntxtbits = 4;
  config.Nuwbits = 5*config.bps;
  config.ftwindow_width = 32;
  config.timing_mx_thresh  = 0.35;
  config.bad_uw_errors = 3;
  config.amp_scale = 245E3;
  config.amp_est_mode = 0;
  config.EsNo_est_all_symbols = 1;
  config.EsNodB = 3;
  config.state_machine = "voice1";
  config.edge_pilots = 1;
  config.clip_gain1 = 2.5;
  config.clip_gain2 = 0.8;
  config.foff_limiter = 0;
  config.txbpf_width_Hz = 2000;
  config.data_mode = "";

  if strcmp(mode,"700D") ||  strcmp(mode,"700d")
    % defaults above
  elseif strcmp(mode,"700E") ||  strcmp(mode,"700e")
    Ts = 0.014; Tcp=0.006; Nc = 21; Ns=4;
    config.edge_pilots = 0; config.state_machine = "voice2";
    config.Nuwbits = 12; config.bad_uw_errors = 3; config.Ntxtbits = 2;
    config.amp_est_mode = 1; config.ftwindow_width = 80;
    config.amp_scale = 155E3; config.clip_gain1 = 3; config.clip_gain2 = 0.8;
    config.foff_limiter = 1;
  elseif strcmp(mode,"2020")
    Ts = 0.0205; Nc = 31;
    config.amp_scale = 167E3; config.clip_gain1 = 2.5; config.clip_gain2 = 0.8;
  elseif strcmp(mode,"2020B")
    Ts = 0.014; Tcp = 0.004; Nc = 29; Ns=5;
    config.Ntxtbits = 4; config.Nuwbits = 8*2; config.bad_uw_errors = 5;
    config.amp_scale = 130E3; config.clip_gain1 = 2.5; config.clip_gain2 = 0.8;
    config.edge_pilots = 0; config.state_machine = "voice2";
    config.foff_limiter = 1; config.ftwindow_width = 64;
    config.txbpf_width_Hz = 2200;
  elseif strcmp(mode,"qam16c1")
    Ns=5; config.Np=5; Tcp = 0.004; Ts = 0.016; Nc = 33; config.data_mode = "streaming";
    config.bps=4; config.Ntxtbits = 0; config.Nuwbits = 15*4; config.bad_uw_errors = 5;
    config.state_machine = "data";
    config.ftwindow_width = 32; config.amp_scale = 132E3;
    config.EsNo_est_all_symbols = 0; config.amp_est_mode = 1; config.EsNodB = 10;
  elseif strcmp(mode,"qam16c2")
    Ns=5; config.Np=31; Tcp = 0.004; Ts = 0.016; Nc = 33; config.data_mode = "streaming";
    config.bps=4; config.Ntxtbits = 0; config.Nuwbits = 42*4; config.bad_uw_errors = 15;
    config.ftwindow_width = 80; config.amp_scale = 135E3; config.state_machine = "data";
    config.EsNo_est_all_symbols = 0; config.amp_est_mode = 1; config.EsNodB = 10;
    config.tx_uw = zeros(1,config.Nuwbits = 42*4);
    config.tx_uw(1:24) = [1 1 0 0  1 0 1 0  1 1 1 1  0 0 0 0  1 1 1 1  0 0 0 0];
    config.tx_uw(end-24+1:end) = [1 1 0 0  1 0 1 0  1 1 1 1  0 0 0 0  1 1 1 1  0 0 0 0];
  elseif strcmp(mode,"datac0")
    Ns=5; config.Np=4; Tcp = 0.006; Ts = 0.016; Nc = 9; config.data_mode = "streaming";
    config.Ntxtbits = 0; config.Nuwbits = 32; config.bad_uw_errors = 9;
    config.state_machine = "data";
    config.ftwindow_width = 80; config.amp_est_mode = 1; config.EsNodB = 3;
    config.edge_pilots = 0; config.timing_mx_thresh = 0.08;
    config.tx_uw = zeros(1,config.Nuwbits);
    config.tx_uw(1:16) = [1 1 0 0  1 0 1 0  1 1 1 1  0 0 0 0];
    config.amp_scale = 300E3; config.clip_gain1 = 2.2; config.clip_gain2 = 0.85;
  elseif strcmp(mode,"datac5")
    Ns=5; config.Np=58; Tcp = 0.004; Ts = 0.016; Nc = 35; config.data_mode = "streaming";
    config.Ntxtbits = 0; config.Nuwbits = 40; config.bad_uw_errors = 14;
    config.state_machine = "data";
    config.ftwindow_width = 80; config.amp_est_mode = 1; config.EsNodB = 3;
    config.amp_scale = 145E3; config.clip_gain1 = 2.7; config.clip_gain2 = 0.8;
    config.edge_pilots = 0; config.timing_mx_thresh = 0.10;
    config.tx_uw = zeros(1,config.Nuwbits);
    config.tx_uw(1:16) = [1 1 0 0  1 0 1 0  1 1 1 1  0 0 0 0];
  elseif strcmp(mode,"datac1")
    Ns=5; config.Np=38; Tcp = 0.006; Ts = 0.016; Nc = 27; config.data_mode = "streaming";
    config.Ntxtbits = 0; config.Nuwbits = 16; config.bad_uw_errors = 6;
    config.state_machine = "data";
    config.ftwindow_width = 80; config.amp_est_mode = 1; config.EsNodB = 3;
    % clipper/compression adjustment:
    % 1. With clipper off increase amp_scale until peak just hit 16384
    % 2. With clipper on increase clip_gain1 until about 30% clipped
    % 3. BPF will drop level beneath 16384, adjust clip_gain2 to just hit 16384 peak again
    % 4. Clipped/unclipped operating point for same PER should be about 1dB apart
    config.amp_scale = 145E3; config.clip_gain1 = 2.7; config.clip_gain2 = 0.8;
    config.edge_pilots = 0; config.timing_mx_thresh = 0.10;
    config.tx_uw = [1 1 0 0  1 0 1 0  1 1 1 1  0 0 0 0];
  elseif strcmp(mode,"datac3")
    Ns=5; config.Np=29; Tcp = 0.006; Ts = 0.016; Nc = 9; config.data_mode = "streaming";
    config.edge_pilots = 0;
    config.Ntxtbits = 0; config.Nuwbits = 40; config.bad_uw_errors = 10;
    config.ftwindow_width = 80; config.timing_mx_thresh = 0.10;
    config.tx_uw = zeros(1,config.Nuwbits);
    config.tx_uw(1:24) = [1 1 0 0  1 0 1 0  1 1 1 1  0 0 0 0  1 1 1 1  0 0 0 0];
    config.tx_uw(end-24+1:end) = [1 1 0 0  1 0 1 0  1 1 1 1  0 0 0 0  1 1 1 1  0 0 0 0];
    config.amp_est_mode = 1; config.EsNodB = 3;
    config.state_machine = "data"; 
    config.amp_scale = 300E3; config.clip_gain1 = 2.2; config.clip_gain2 = 0.8;
  elseif strcmp(mode,"datac4")
    Ns=5; config.Np=47; Tcp = 0.006; Ts = 0.016; Nc = 4; config.data_mode = "streaming";
    config.edge_pilots = 0;
    config.Ntxtbits = 0; config.Nuwbits = 32; config.bad_uw_errors = 12;
    config.ftwindow_width = 80; config.timing_mx_thresh = 0.5;
    config.tx_uw = zeros(1,config.Nuwbits);
    config.tx_uw(1:24) = [1 1 0 0  1 0 1 0  1 1 1 1  0 0 0 0  1 1 1 1  0 0 0 0];
    config.tx_uw(end-24+1:end) = [1 1 0 0  1 0 1 0  1 1 1 1  0 0 0 0  1 1 1 1  0 0 0 0];
    config.amp_est_mode = 1; config.EsNodB = 3;
    config.state_machine = "data";
    config.amp_scale = 2*300E3; config.clip_gain1 = 1.2; config.clip_gain2 = 1.0;
    config.txbpf_width_Hz = 400;
 elseif strcmp(mode,"datac13")
    Ns=5; config.Np=18; Tcp = 0.006; Ts = 0.016; Nc = 3; config.data_mode = "streaming";
    config.edge_pilots = 0;
    config.Ntxtbits = 0; config.Nuwbits = 48; config.bad_uw_errors = 18;
    config.ftwindow_width = 80; config.timing_mx_thresh = 0.45;
    config.tx_uw = zeros(1,config.Nuwbits);
    config.tx_uw(1:24) = [1 1 0 0  1 0 1 0  1 1 1 1  0 0 0 0  1 1 1 1  0 0 0 0];
    config.tx_uw(end-24+1:end) = [1 1 0 0  1 0 1 0  1 1 1 1  0 0 0 0  1 1 1 1  0 0 0 0];
    config.amp_est_mode = 1; config.EsNodB = 3;
    config.state_machine = "data"; 
    config.amp_scale = 2.5*300E3; config.clip_gain1 = 1.2; config.clip_gain2 = 1.0;
    config.txbpf_width_Hz = 400;
  elseif strcmp(mode,"1")
    Ns=5; config.Np=10; Tcp=0; Tframe = 0.1; Ts = Tframe/Ns; Nc = 1;
  else
    % try to parse mode string for user defined mode
    vec = sscanf(mode, "Ts=%f Nc=%d Ncp=%f");
    Ts=vec(1); Nc=vec(2); Ncp=vec(3);
  end
  Rs=1/Ts;
  config.Rs = Rs; config.Tcp = Tcp; config.Ns = Ns; config.Nc = Nc;
  if !isfield(config,"tx_uw") 
    config.tx_uw = zeros(1,config.Nuwbits); 
  end  
end

% ------------------------------------------------------------------------------
% codec_to_frame_packing - Set up a bunch of constants to support modem frame
%                          construction from LDPC codewords and codec source bits
% ------------------------------------------------------------------------------

function [code_param Nbitspercodecframe Ncodecframespermodemframe] = codec_to_frame_packing(states, mode)
  ofdm_load_const;
  mod_order = 4; bps = 2; modulation = 'QPSK'; mapping = 'gray';

  init_cml();
  if strcmp(mode, "700D")
    load HRA_112_112.txt
    code_param = ldpc_init_user(HRA_112_112, modulation, mod_order, mapping);
    assert(Nbitsperframe == (code_param.coded_bits_per_frame + Nuwbits + Ntxtbits));
    % unused for this mode
    Nbitspercodecframe = Ncodecframespermodemframe = 0;
  end
  if strcmp(mode, "700E")
    load HRA_56_56.txt
    code_param = ldpc_init_user(HRA_56_56, modulation, mod_order, mapping);
    assert(Nbitsperframe == (code_param.coded_bits_per_frame + Nuwbits + Ntxtbits));
    % unused for this mode
    Nbitspercodecframe = Ncodecframespermodemframe = 0;
  end
  if strcmp(mode, "2020")
    load HRA_504_396.txt
    code_param = ldpc_init_user(HRA_504_396, modulation, mod_order, mapping);
    code_param.data_bits_per_frame = 312;
    code_param.coded_bits_per_frame = code_param.data_bits_per_frame + code_param.ldpc_parity_bits_per_frame;
    code_param.coded_syms_per_frame = code_param.coded_bits_per_frame/code_param.bits_per_symbol;
    printf("2020 mode\n");
    printf("ldpc_data_bits_per_frame = %d\n", code_param.ldpc_data_bits_per_frame);
    printf("ldpc_coded_bits_per_frame  = %d\n", code_param.ldpc_coded_bits_per_frame);
    printf("ldpc_parity_bits_per_frame  = %d\n", code_param.ldpc_parity_bits_per_frame);
    printf("data_bits_per_frame = %d\n", code_param.data_bits_per_frame);
    printf("coded_bits_per_frame  = %d\n", code_param.coded_bits_per_frame);
    printf("coded_syms_per_frame  = %d\n", code_param.coded_syms_per_frame);
    printf("ofdm_bits_per_frame  = %d\n", Nbitsperframe);
    Nbitspercodecframe = 52; Ncodecframespermodemframe = 6;
    printf("  Nuwbits: %d  Ntxtbits: %d\n", Nuwbits, Ntxtbits);
    Nparity = code_param.ldpc_parity_bits_per_frame;
    totalbitsperframe = code_param.data_bits_per_frame + Nparity + Nuwbits + Ntxtbits;
    printf("Total bits per frame: %d\n", totalbitsperframe);
    assert(totalbitsperframe == Nbitsperframe);
  end
  if strcmp(mode, "qam16c1")
      load H2064_516_sparse.mat
      code_param = ldpc_init_user(HRA, modulation='QAM', mod_order=16, mapping="", reshape(states.qam16,1,16));
  end
  if strcmp(mode, "qam16c2")
      framesize = 16200; rate = 0.6;
      code_param = ldpc_init_builtin("dvbs2", rate, framesize, modulation='QAM', mod_order=16, mapping="", reshape(states.qam16,1,16));
  end
  if strcmp(mode, "datac5")
      framesize = 16200; rate = 0.6;
      code_param = ldpc_init_builtin("dvbs2", rate, framesize, modulation='QPSK', mod_order=4, mapping="");
  end
  if strcmp(mode, "datac0") || strcmp(mode, "datac13")
    load H_128_256_5.mat
    code_param = ldpc_init_user(H, modulation, mod_order, mapping);
  end
  if strcmp(mode, "datac1")
    load H_4096_8192_3d.mat
    code_param = ldpc_init_user(HRA, modulation, mod_order, mapping);
  end
  if strcmp(mode, "datac3")
    load H_1024_2048_4f.mat
    code_param = ldpc_init_user(H, modulation, mod_order, mapping);
  end
  if strcmp(mode, "datac4")
    load H_1024_2048_4f
    code_param = ldpc_init_user(H, modulation, mod_order, mapping);
    code_param.data_bits_per_frame = 448;
    code_param.coded_bits_per_frame = code_param.data_bits_per_frame + code_param.ldpc_parity_bits_per_frame;
    code_param.coded_syms_per_frame = code_param.coded_bits_per_frame/code_param.bits_per_symbol;
  end
  if strcmp(mode, "datac13")
    load H_256_512_4.mat
    code_param = ldpc_init_user(H, modulation, mod_order, mapping);
    code_param.data_bits_per_frame = 128;
    code_param.coded_bits_per_frame = code_param.data_bits_per_frame + code_param.ldpc_parity_bits_per_frame;
    code_param.coded_syms_per_frame = code_param.coded_bits_per_frame/code_param.bits_per_symbol;
  end
  if strcmp(mode, "datac0") || strcmp(mode, "datac1") || strcmp(mode, "datac3") ...
     || strcmp(mode, "datac4") || strcmp(mode, "qam16c1") ...
     || strcmp(mode, "qam16c2") || strcmp(mode, "datac5") || strcmp(mode, "datac13")
    printf("ldpc_data_bits_per_frame = %d\n", code_param.ldpc_data_bits_per_frame);
    printf("ldpc_coded_bits_per_frame  = %d\n", code_param.ldpc_coded_bits_per_frame);
    printf("ldpc_parity_bits_per_frame  = %d\n", code_param.ldpc_parity_bits_per_frame);
    printf("Nbitsperpacket  = %d\n", Nbitsperpacket);
    Nparity = code_param.ldpc_parity_bits_per_frame;
    totalbitsperframe = code_param.data_bits_per_frame + Nparity + Nuwbits + Ntxtbits;
    printf("totalbitsperframe = %d\n", totalbitsperframe);
    assert(totalbitsperframe == Nbitsperpacket);
    Nbitspercodecframe = Ncodecframespermodemframe = -1;
  end
endfunction


