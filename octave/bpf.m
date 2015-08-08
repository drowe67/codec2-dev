% bpf.m
% David Rowe April 2015
%
% Design 400-2600 Hz BPF and save coeffs 

1;

function write_c_array(filename, arrayname, vec)

  m = length(vec);

  f=fopen(filename,"wt");
  fprintf(f,"#define %s_N %d\n\n", toupper(arrayname), m);
  fprintf(f,"float %s[]={\n", arrayname);
  for r=1:m
    if r < m
      fprintf(f, "  %f,\n", vec(r));
    else
      fprintf(f, "  %f\n};", vec(r));
    end
  end

  fclose(f);
endfunction

b=firls(100,[0 250 300 2600 2700 4000]/4000,[0.01 0.01 1 1 0.01 0.01]);
freqz(b)
write_c_array("../src/bpfb.h", "bpfb", b)

% C header file of noise samples so C version gives extacly the same results

