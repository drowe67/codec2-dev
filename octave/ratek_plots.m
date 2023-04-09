% ratek1_plots.m
%
% David Rowe 2022
%
% Plots for Rate K report.
%

function ratek_plots(epslatex_path)
  more off;

  newamp_700c;
  Fs = 8000;  K = 20;

  figure(1); clf; 
  f = 0:4000; plot(f, ftomel(f));
  grid;
  set(gca, 'FontSize', 16);
  xlabel('Frequency f (Hz)'); ylabel('mel(f)');
  print("ratek_mel_fhz","-dpng","-S500,500");

  figure(1); clf; 
  k = 1:K; plot(k, warp(k,K));
  grid;
  set(gca, 'FontSize', 16);
  xlabel('k'); ylabel('f (Hz)');
  print("warp_fhz_k","-dpng","-S500,500");

  figure(1); clf; 
  Fs=8000; f0=200; L=floor(Fs/(2*f0)); Nb=10;
  hold on;
  for m=1:L
    plot((1:L)*f0,generate_filter(m,f0,L,Nb),'o-');
  end
  hold off; grid;
  set(gca, 'FontSize', 16);
  xlabel('f (Hz)'); ylabel('h(m=f/F0)');
  print("filters_h","-dpng","-S500,500");

  % show effect of Hilbert Compressor and energy limiting
  figure(1); clf;
  if length(epslatex_path)
    textfontsize = get(0,"defaulttextfontsize");
    linewidth = get(0,"defaultlinelinewidth");
    set(0, "defaulttextfontsize", 10);
    set(0, "defaultaxesfontsize", 10);
    set(0, "defaultlinelinewidth", 0.5);
  end
  s1=load_raw("../build_linux/230331_mean_energy/forig_5_vq1_12_dec3.wav");
  s2=load_raw("../build_linux/230331_mean_energy/forig_6_vq1_12_dec3_q4.wav");
  s3=load_raw("../build_linux/230331_mean_energy/forig_7_vq1_12_dec3_q4_hc.wav");
  subplot(311); plot(s1); axis([100 length(s1) -3.2E4 3.2E4]); grid;
  subplot(312); plot(s2); axis([100 length(s2) -3.2E4 3.2E4]); grid;
  subplot(313); plot(s3); axis([100 length(s3) -3.2E4 3.2E4]); grid;
  if length(epslatex_path)
    old_dir=cd(epslatex_path);
    fn = "hilbert_clipper.tex"
    print(fn,"-depslatex","-S400,400");
    printf("printing... %s%s\n", epslatex_path,fn);
    cd(old_dir);
    set(0, "defaulttextfontsize", textfontsize);
    set(0, "defaultaxesfontsize", textfontsize);
    set(0, "defaultlinelinewidth", linewidth);
  end
endfunction
