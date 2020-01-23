function assert_approx_eq(actual, expected, epsilon)

if (abs(actual - expected) > epsilon)
   error('Values are not equal! Expected: %d, actual: %d, diff: %d',...
       expected, actual, abs(actual - expected));    
end

end