function ke_met_final = ke_met_solve2C(parameters, ode, y0, ke)
tol = .01; %this is the tolerance of the solver. It is a fraction of the literatire half life, thus a tolerance of .01 is 1% of the literature half life.
%The half life returned by the guessed Ke value must be within the
%tolerance to be considered correct


%the half-life of the metabolite (hr)
t_half_met = parameters(10);


T_half = t_half_met/2; %need to define some arbitrary experimental T_half for loop to work
    low_guess = 0;
    high_guess = 10*log(2)/t_half_met;
    %loop runs while the experimental T_Half it outside the tolerance
while T_half < t_half_met - tol * t_half_met || T_half > t_half_met + tol * t_half_met
    ke_met = (low_guess + high_guess)/2; %Ke is guessed halfway between the bounds
    if ke == 0 || ke_met == 10*log(2)/t_half_met
        error('Bounds of high and low guess for ke_met_solve2C are not wide enough. Go in and change it.')
    elseif abs(high_guess-low_guess) < 0.0001
         ke_met_final = ke_met;
         break
    end
    build_ode = @(t, x) ode(ke_met, t, x);
    
    [T,Y] = ode45(build_ode, [0 1000], y0);
    
    logC1 = log(Y(:,3));
    j = 1;
    while T(j) <  4*t_half_met
        j = j + 1;
    end %This searches for a timepoint deep into the profile to avoid distribution effects
    j2 = 1;
    while T(j2) < 5*t_half_met
        j2 = j2 + 1;
    end
    q = polyfit(T(j:j2),logC1(j:j2),1); %fits a line to the log plot slope
    T_half = -log(2)/q(1); %the slope of the line is -ln(2) / t_half_met, solve for guess t_half_met
    
    
    if T_half > t_half_met - tol * t_half_met %adjust the upper and lower bounds based on if the guess was too high or too low
        low_guess = ke_met;
    elseif T_half < t_half_met + tol * t_half_met
        high_guess = ke_met;
    end
end
ke_met_final = ke_met;
end
