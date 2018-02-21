function [estimates, model] = fitCurve(data, start_point)
% Call fminsearch with a random starting point.
if ~exist('start_point') start_point = [0.5 0.001 0.7  ];  end %rand(1, 3);
xdata=data(:,1);ydata=data(:,2);
model = @expfun;
options=optimset('display','off','Funvalcheck','off','maxFunevals',10000,'maxiter',10000);
[estimates,fval,exitflag,output] = fminsearch(model, start_point,options);
% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A * exp(-lambda * xdata) - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    function [sse, FittedCurve] = expfun(params)
        k1 = params(1);
        k2 = params(2);
        k3 = params(3);
        FittedCurve = (k1+k2*xdata).*(1-exp(-k3*xdata));
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end
end