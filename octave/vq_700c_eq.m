% vq_700c_plots
% David Rowe May 2019

function vq_700c_plots(fn_array)
  nb_features = 41
  K = 20
  figure(1); clf; hold on; axis([1 20 -20 40]); title('Max Hold');
  figure(2); clf; hold on; axis([1 20 -20 30]); title('Average'); 
  for i=1:length(fn_array)
    [dir name ext] = fileparts(fn_array{i});
    if strcmp(ext, ".f32")
      % f32 feature file
      fn = sprintf("../script/%s_feat%s", name, ext)
      feat = load_f32(fn , nb_features);
      bands = feat(:,2:K+1);
    else
      % text file (e.g. existing VQ)
      bands = load(fn_array{i});
    end
    figure(1); plot(max(bands),'linewidth', 5);
    figure(2); plot(mean(bands),'linewidth', 5);
  end
  figure(1); legend(fn_array);
  figure(2); legend(fn_array);
endfunction
