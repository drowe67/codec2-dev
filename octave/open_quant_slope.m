function b = open_quant_slope(b)
  b(2) = quantise([-1 -0.5 0.5 1], b(2));
  b(1) = quantise([0.5 1.0], b(1));
end
