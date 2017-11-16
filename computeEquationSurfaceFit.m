function [ fitresult, gof, zTeSim, rmseTe ] = computeEquationSurfaceFit( x, y, z, xTe, yTe, zTe, equation, flagShow )
%Author: ylonge.
%Function: compute fit parameter and show curve with given equation. Only one independent factor is considered.
%-Input:
%   --x,y: independent and dependent data.
%	--equation: string that describes the fit equation.
%-OutPut:
%	--fitresult: the parameters of the equation.
%	--gof: performance evaluation of the fit.

%% Fit.
[xData, yData, zData] = prepareSurfaceData( x, y, z );

% Set up fittype and options.
ft = fittype( equation );

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft );

% Compute test data.
zTeSim = fitresult(xTe,yTe);
rmseTe = (mean((zTe - zTeSim).^2))^0.5;

if flagShow
	% Plot fit with data.
	figure( 'Name', 'fit surface' );
% 	h = plot( fitresult, [xData, yData], zData );
%     hold on;    
%     plot3( xTe, yTe, zTe, 'rd', 'MarkerFaceColor','r');
    
    h = plot( fitresult, [xTe, yTe], zTe );
    hold on
    plot3( xData, yData, zData, 'rd', 'MarkerFaceColor','r');
    
	legend( h, 'fit surface', 'z vs. x, y', 'Location', 'NorthEast' );
    % Label axes
    xlabel x
    ylabel y
    zlabel z
    grid on
    view( -0.7, 38.0 );
    axis([min([xTe;xData]) max([xTe;xData]) min([yTe;yData]) max([yTe;yData]) -15 5]);
end
end