function [J, grad] = linearRegCostFunction(X, y, theta, lambda)
%LINEARREGCOSTFUNCTION Compute cost and gradient for regularized linear 
%regression with multiple variables
%   [J, grad] = LINEARREGCOSTFUNCTION(X, y, theta, lambda) computes the 
%   cost of using theta as the parameter for linear regression to fit the 
%   data points in X and y. Returns the cost in J and the gradient in grad

% Initialize some useful values
m = length(y); % number of training examples
h = X * theta;
non_reg_error = sum((h - y).^2)/m/2;
reg_error = sum(theta(2:end).^2)  * lambda /2/m;

% You need to return the following variables correctly 
J = non_reg_error + reg_error;

grad_non_reg = (h - y)' * X / m;
grad_reg = [0; theta(2:end)(:) * lambda / m];

grad = grad_non_reg(:) + grad_reg(:);

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost and gradient of regularized linear 
%               regression for a particular choice of theta.
%
%               You should set J to the cost and grad to the gradient.
%












% =========================================================================

grad = grad(:);

end
