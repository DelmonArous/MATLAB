function [output_arg] = f_1(t)
    output_arg = (pi^2 - t.^2).*exp(t./pi);
end
