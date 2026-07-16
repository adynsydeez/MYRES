function [D, err] = fit_voxel_mono(Svec, bvec19, Dmin, Dmax)
    if Svec(1) <= 0
        D = 0; err = 0; return;
    end
    
    D_custom_guess = 1e-3;

    options = optimoptions('lsqcurvefit', 'Display', 'off', 'FunctionTolerance', 1e-5);
    model = @(p, b) Svec(1) * exp(-b * p(1));
    
    [p, resnorm] = lsqcurvefit(model, D_custom_guess, bvec19, Svec, Dmin, Dmax, options);
    
    D = p(1);
    err = sqrt(resnorm / length(bvec19)); % Percentage error for clim([0,20])
end