function At = ctranspose(A)

A.adjoint = xor(A.adjoint,true);
At = A;
