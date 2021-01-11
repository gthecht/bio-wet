function x_dot = get_dot(a, b, x)
    %get_dot Returns the dot of the letter that has alpha and beta
	x_dot = a .* (1 - x) - b .* x;
end

