function D = difference_operator(n, order)
%DIFFERENCE_OPERATOR Sparse finite-difference operator of a given order.
%   D = DIFFERENCE_OPERATOR(N, ORDER) is the (N - ORDER)-by-N sparse matrix
%   whose rows apply the ORDER-th forward difference to a vector of length
%   N: first order [-1 1], second order [1 -2 1], and so on (signs
%   alternate from (-1)^ORDER). ORDER 0 gives the identity. It replaces
%   get_l of Hansen's Regularization Tools, with the same entries.

if order < 0
    error("difference_operator:Order", "ORDER must be nonnegative.")
end
if order == 0
    D = speye(n);
    return
end
% Coefficients of (x - 1)^order: the order-th forward difference
coefficients = [-1 1];
for k = 2:order
    coefficients = [0 coefficients] - [coefficients 0];
end
n_row = n - order;
D = sparse(n_row, n);
for k = 1:order + 1
    D = D + sparse(1:n_row, (1:n_row) + k - 1, coefficients(k), n_row, n);
end
end
