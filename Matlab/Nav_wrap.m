% Wraps specified 'value' around specified 'bound'
% 'bound' should be grater than zero

function [ result ] = Nav_wrap(value, bound)

if (bound <= 0.0)
    error('Specified "bound" value should be greater than zero');
end

vl = abs(value);
sign_ = sign(value);

while (vl > bound)
    vl = vl - bound;
end

result = vl * sign_;

end