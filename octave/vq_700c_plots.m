% vq_700c_plots
% David Rowe May 2019

function vq_700c_plots(fn_f32)
  nb_features = 41
  K = 20
  figure(1); clf; hold on; axis([1 20 -20 40]); title('Max Hold');
  figure(2); clf; hold on; axis([1 20 -20 30]); title('Average'); 
  for i=1:length(fn_f32)
    fn = sprintf("../script/%s_feat.f32", fn_f32{i});
    feat = load_f32(fn , nb_features);
    bands = feat(:,2:K+1);
    figure(1); plot(max(bands));
    figure(2); plot(mean(bands));
  end
  figure(1); legend(fn_f32);
  figure(2); legend(fn_f32);
endfunction
