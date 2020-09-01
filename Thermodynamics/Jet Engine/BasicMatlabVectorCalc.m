% Basic matlab knowledge on vector 
a = [1 2 3];b=[4 5 6];
A = a'*b    % A is a matrix
A = a*b'    % A is a scalar
A = a*b     % gives an error message 
            % Error using  * 
            % Inner matrix dimensions must agree.
