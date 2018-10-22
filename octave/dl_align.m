% dl_align.m
% David Rowe Oct 2018
%
% visual test of alignment rqd for deep learning experients

s16 = load_raw("~/Downloads/16k-LP7/CA/c01_01.raw");
s8 = load_raw("../build_linux/src/centre.raw");

% time axis

t16 = 0:1/16000:length(s16)/16000;
t8 = 0:1/8000:length(s8)/8000;

st_t = 1.12; en_t = 1.16;
d_16 = floor(47/2 + 0.015*16000)
printf("delay 16kHz by %d samples\n", d_16)
st_16 = st_t*16000; en_16 = en_t*16000;
st_8 = st_t*8000; en_8 = en_t*8000;

figure(1); clf;
plot(t16(st_16:en_16), s16(st_16-d_16:en_16-d_16));
hold on; plot(t8(st_8:en_8), s8(st_8:en_8),'g');
axis([st_t en_t -3000 3000])
figure(2); clf;
subplot(211); plot(s16);
subplot(212); plot(s8);
