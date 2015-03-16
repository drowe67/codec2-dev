% lspvar.m
% David Rowe 8 Feb 2015
%
% Experiments in reducing LSP variance in order to make quantisation easier

lsp;
close all;
graphics_toolkit ("gnuplot");

function v = det_var(s)
  p = 10;
  ak=lpcauto(s,p,[160 80]);
  [rows cols] = size(ak);
  w = zeros(rows,p);
  for r=1:rows
    w(r,:) = atolsp(ak(r,:));
  end
  v = sum(var(w));
  figure;
  plot(w)
end

s=load_raw("../raw/cq_ref.raw");
det_var(s)

if 1
 Fs = 8000;
 [b a] = cheby2(6,40,[300 2600]/(Fs/2));
 s1 = filter(b,a,s);
 det_var(s1)
end

% + Equalise to remove spectral slope, this reduces variance of spectral 
%   envelope
% + Can probably get rid of energy for many vectors too, except flat ones?
% + Or have some sort of correlation of gain with peakiness of vectors?
% + Once flattened try slicing out low energy harmonics, like post filter.  How 
%   does this sound?  Any artefacts of harmonics coming and going near 
%   thresholds?
% + BPF as well, it's all abt minimal information.
% + What's left?  Look at it.  How can it be represented?  Look at tol distort,
%   what quantiser will deliver that?  Just position of clipped peaks?
% + fall back HF mode, no energy, voicing, just spectral information.  Diveristy
%   compared to FEC?  Test with simulation of ideal code.  Aux carrier with 
%   extra information to get 1300-ish quality.  Coherent, pilot assisted demod.
%   Low PAPR as extra carriers at lower power level.  Essentially for free in
%   terms of TX power.
% + What is a viable plan to get this going in a mon month or so?
%   + Octave equlaisation of input files
%   + take model files, simulate quantisation, back to C for synth?
%   + Coherent FDM modem port to C






