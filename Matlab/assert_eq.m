function assert_eq(actual, expected, message)

if (actual ~= expected)
   error(message);
end

end