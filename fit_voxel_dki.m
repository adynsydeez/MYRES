function [D, K, err] = fit_voxel_dki(Svec, bvec, D_guess, K_guess, Dmin, Dmax, Kmin, Kmax)
    if Svec(1) <= 0
        D = 0; K = 0; err = 0; return;
    end
    
    options = optimoptions('lsqcurvefit', 'Display', 'off', 'FunctionTolerance', 1e-5);
    
    model = @(p, b) Svec(1) * exp(-b * p(1) + (1/6) * (b.^2) * (p(1)^2) * p(2));
    
    [p, resnorm] = lsqcurvefit(model, [D_guess, K_guess], bvec, Svec, [Dmin, Kmin], [Dmax, Kmax], options);
    
    D = p(1);
    K = p(2);
    
    err = sqrt(resnorm / length(bvec));
end