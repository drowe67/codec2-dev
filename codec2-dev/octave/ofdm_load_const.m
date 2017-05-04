% make like C #define for ofdm modem

Fs = states.Fs;
bps = states.bps;
Rs = states.Rs;
Tcp = states.Tcp;
Ns = states.Ns;
Nc = states.Nc;
M = states.M;
Ncp = states.Ncp;
bps = states.bps;
Nbitsperframe = states.Nbitsperframe;
Nrowsperframe = states.Nrowsperframe;
Nsamperframe = states.Nsamperframe;

W = states.W;
w = states.w;

timing_est = states.timing_est;
sample_point = states.sample_point;
ftwindow_width = states.ftwindow_width;

Nrxbuf = states.Nrxbuf;
rxbuf = states.rxbuf;

pilots = states.pilots;
rate_fs_pilot_samples = states.rate_fs_pilot_samples;

foff_est_gain = states.foff_est_gain;
foff_est_hz = states.foff_est_hz;

timing_en = states.timing_en;
foff_est_en = states.foff_est_en;
phase_est_en = states.phase_est_en;

rate = states.rate;
ldpc_en = states.ldpc_en;
if ldpc_en
  code_param = states.code_param;
  max_iterations = states.ldpc_max_iterations;
  demod_type = states.ldpc_demod_type;
  decoder_type = states.ldpc_decoder_type;
end

verbose = states.verbose;

