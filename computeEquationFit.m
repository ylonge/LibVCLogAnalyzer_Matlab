function [ fitresult, gof ] = computeEquationFit( x, y, equation, startPoint, flagShow )
%Author: ylonge.
%Function: compute fit parameter and show curve with given equation. Only one independent factor is considered.
%-Input:
%   --x,y: independent and dependent data.
%	--equation: string that describes the fit equation.
%-OutPut:
%	--fitresult: the parameters of the equation.
%	--gof: performance evaluation of the fit.

%% Fit.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( equation, 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = startPoint;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if flagShow
	% Plot fit with data.
	figure( 'Name', 'fit curve' );
	h = plot( fitresult, xData, yData );
	legend( h, 'y vs. x', 'fit curve', 'Location', 'NorthEast' );
	% Label axes
	xlabel( 'x' );
	ylabel( 'y' );
	grid on
end
end