function ka_final = ka_solve2C(parameters, ode, y0)
%Ka is the same for all oral systems
tol = .1;
low_guess_ka = 0; %define some high and low bounds for guessing
high_guess_ka = 100;
%The time where max concentration occurs
t_max = parameters(5);

T_max = t_max/2;
while T_max < t_max - tol * t_max || T_max > t_max + tol * t_max
    ka = (low_guess_ka + high_guess_ka)/2;
    if ka == 0 || ka == 100
        error('Ka solver hit the bounds. Go to extended1C_noMet.m to change low and high guess.')
    end
    build_ode = @(t, x) ode(ka, t, x);
    [T, Y] = ode45(build_ode, [0 1000], y0);
    [argmax, maxVal] = max(Y(:,2));
    T_max = T(maxVal);
    if T_max < t_max - tol * t_max
        high_guess_ka = ka;
    elseif T_max > t_max + tol * t_max
        low_guess_ka = ka;
    end
end
ka_final = ka;
end