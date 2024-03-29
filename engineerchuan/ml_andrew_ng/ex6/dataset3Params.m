function [C, sigma] = dataset3Params(X, y, Xval, yval)
%DATASET3PARAMS returns your choice of C and sigma for Part 3 of the exercise
%where you select the optimal (C, sigma) learning parameters to use for SVM
%with RBF kernel
%   [C, sigma] = DATASET3PARAMS(X, y, Xval, yval) returns your choice of C and 
%   sigma. You should complete this function to return the optimal C and 
%   sigma based on a cross-validation set.
%

% You need to return the following variables correctly.
C = 1;
sigma = 0.3;

% ====================== YOUR CODE HERE ======================
% Instructions: Fill in this function to return the optimal C and sigma
%               learning parameters found using the cross validation set.
%               You can use svmPredict to predict the labels on the cross
%               validation set. For example, 
%                   predictions = svmPredict(model, Xval);
%               will return the predictions on the cross validation set.
%
%  Note: You can compute the prediction error using 
%        mean(double(predictions ~= yval))
%
params = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30];
%params = [0.01, 0.03, 0.1];
values = zeros(length(params) * length(params), 3);
for i_C = 1: length(params)
  for i_sigma = 1:length(params)
    C = params(i_C);
    sigma = params(i_sigma);
    model= svmTrain(X, y, C, @(x1, x2) gaussianKernel(x1, x2, sigma));
    pred_values = svmPredict(model, Xval);
    error = mean(double(pred_values ~= yval));
    values((i_C-1) * length(params) + i_sigma, :) = [error, C, sigma];
  endfor
endfor
values;
% find the max
[min_val_error, min_val_error_index] = min(values(:, 1));
C = values(min_val_error_index, 2)
sigma = values(min_val_error_index, 3)


end
