% tmic_eq.m
% David Rowe Feb 2023
%
% Another attempt at microphone equaliser

#{
  cd ~/codec2/build_linux
  ./src/c2sim --hpf ../raw/big_dog.raw --modelout big_dog_model.bin
#}

function tmic_eq(samname, vq_stage1_f32, frame_limit)
  K = 20;
  newamp_700c;
  
  model_name = strcat(samname,"_model.bin");
  model = load_codec2_model(model_name);
  [frames tmp] = size(model);
  if nargin == 3
    frames = frame_limit
    model = model(1:frames,:);
  end
  rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K);
  vq_stage1 = load_f32(vq_stage1_f32,K);
  vq(:,:,1)= vq_stage1;
  printf("Building input B vectors ...\n");
  B = ratek3_batch_tool(samname,'K',K);
  B = B(1:frames,:);
  
  printf("Training mic equaliser...\n");
  f = zeros(1,K); f_log = zeros(frames,K);
  for i=1:frames
    f = eq_cand_b(f, B(i,:), vq, delta=0.01);
    f_log(i,:) = f;
    printf("%d/%d %3.0f%%\r", i,frames, (i/frames)*100)
  end
  printf("\n");
  
  figure(1); clf; 
  plot(rate_K_sample_freqs_kHz*1000, f, 'b;f cand B EQ;');
  q = mean(B-mean(B,2)) - mean(vq_stage1); q = max(q,0);
  hold on; plot(rate_K_sample_freqs_kHz*1000, q, 'r;q newamp1 EQ;'); hold off;
  axis([0 4000 -10 10]); grid('minor');
  figure(2); clf; mesh(f_log);

  % determine VQ distortion without and with equaliser
  printf("Computing distortion...\n");
  target = B - mean(B,2);
  [res target_ ind] = mbest(vq, target, 1);
  Eq_vanilla = mean((target(:) - target_(:)) .^2);
  printf("Eq_vanilla..: %5.2f dB^2\n", Eq_vanilla);
  target = B-q;
  target -= mean(target,2);
  [res target_ ind] = mbest(vq, target, 1);
  Eq_equal_q = mean((target(:) - target_(:)) .^2);
  printf("Eq_newamp1..: %5.2f dB^2\n", Eq_equal_q);
  target = B-f;
  target -= mean(target,2);
  [res target_ ind] = mbest(vq, target, 1);
  Eq_equal_f = mean((target(:) - target_(:)) .^2);
  printf("Eq_candB....: %5.2f dB^2\n", Eq_equal_f);
end
