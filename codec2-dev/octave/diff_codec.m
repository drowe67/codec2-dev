% diff_codec.m
%
% Plots differences between two states in two runs of the codec,
% e.g. x86 and embedded.
%
% Copyright David Rowe 2013
%
% This program is distributed under the terms of the GNU General Public License 
% Version 2

function diff_codec(samname1, samname2, model1_prefix, model2_prefix)
  
  fs1=fopen(samname1,"rb");
  s1=fread(fs1,Inf,"short");
  fs2=fopen(samname2,"rb");
  s2=fread(fs2,Inf,"short");

  st = 1;
  en = length(s1);

  figure(1);
  clf;
  subplot(211);
  l1 = strcat("r;",samname1,";");
  plot(s1(st:en), l1);
  axis([1 en-st min(s1(st:en)) max(s1(st:en))]);
  subplot(212);
  l2 = strcat("r;",samname2,";");
  plot(s2(st:en),l2);
  axis([1 en-st min(s1(st:en)) max(s1(st:en))]);
 
  figure(2)
  plot(s1(st:en)-s2(st:en));
  max(s1(st:en)-s2(st:en));
  
  model_name1 = strcat(model1_prefix,"_model.txt");
  model1 = load(model_name1);
  model_name1q = strcat(model1_prefix,"_qmodel.txt");
  model1q = load(model_name1q);

  model_name2 = strcat(model2_prefix,"_model.txt");
  model2 = load(model_name2);
  model_name2q = strcat(model2_prefix,"_qmodel.txt");
  model2q = load(model_name2q);

  Wo1 = model1(:,1);
  L1 = model1(:,2);
  Am1 = model1(:,3:82);
  Wo1q = model1q(:,1);
  L1q = model1q(:,2);
  Am1q = model1q(:,3:82);

  Wo2 = model2(:,1);
  L2 = model2(:,2);
  Am2 = model2(:,3:82);
  Wo2q = model2q(:,1);
  L2q = model2q(:,2);
  Am2q = model2q(:,3:82);

  figure(3)
  subplot(211)
  plot(Wo1)
  title('Wo1');
  subplot(212)
  plot(Wo1-Wo2)
  figure(4)
  subplot(211)
  plot(Wo1q)
  title('Wo1q');
  subplot(212)
  plot(Wo1q-Wo2q)

  figure(5)
  subplot(211)
  plot(L1)
  title('L1');
  subplot(212)
  plot(L1-L2)
  figure(6)
  subplot(211)
  plot(L1q)
  title('L1q');
  subplot(212)
  plot(L1q-L2q)
  
  figure(7)
  l=length(L1q);
  sm=zeros(1,l);
  for f=1:l
    %printf("f %d L1q %d L2q %d\n",f,L1q(f),L2q(f));  
    sm(f) = sum(10*log10(Am1q(f,1:L1q(f))) - 10*log10(Am2q(f,1:L2q(f))));
  end
  plot(sm)
  title('Am1q - Am2q');

endfunction
