%  Wraps specified 'valuee' around 2pi
function [ result ] = Nav_wrap_2pi(value)
    result = Nav_wrap(value, 2.0 * pi);
end