clear all 
clc
x = linspace(0, 1000, 1001);
%% Just some examples of function definition
% as an algebraic transformation
f = x;
f(1);

% inline function
g = @(x) x;
g(0);

%general function 
function x = h (a)
    x = a; % Calculate the result based on the function f and parameters
end
h(0);