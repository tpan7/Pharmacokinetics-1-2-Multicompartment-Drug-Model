function multiple1C_Met(parameters, route, num_doses)
%Which route of administration
%1 for Intravenous, 2 for Oral, 3 for Intramusclar, 4 for Subcutaneous
%5 for Buccal, 6 for Sublingual, 7 for Lung, 8 for Rectal, 9 for
%transdermal
[temp, indx] = max(route);

%The Volume of distribution
Vd = parameters(1);
%The half life
t_half = parameters(2);
%The bioavailability factor
F = parameters(3);
if indx == 1
    F = 1;
end
%The time where max concentration occurs
t_max = parameters(4);
%The end of the plotting
range = parameters(5);
%the volume of distribution of the metabolite (L/kg)
Vd_met = parameters(6);
%the half-life of the metabolite (hr)
t_half_met = parameters(7);
%the fraction of drug is converted to metabolite
F_met = parameters(8);


%rate of elimination
ke = log(2)/t_half;
%rate of elimination of metabolite
ke_met = log(2)/t_half_met;

dosing_times = [0];
Dose = [];
for i=1:num_doses-1
    ti = input(['At what time did you give dose #' num2str(i+1) ' (hr)? (Assuming dose #1 was at t=0): ']);
    if ti > range
        error('This dose is outside the time you are integrating over.')
    elseif ti <= dosing_times(end)
        error('This dose time is before or at the same time as the previous one.'); 
    end
    dosing_times = [dosing_times ti];
end
for i=1:num_doses
    Dosei = input(['What was dose #' num2str(i) ' (mg)?: ']);
    Dose = [Dose Dosei];
end

%Ka is the same for all oral systems
low_guess_ka = 0; %define some high and low bounds for guessing
high_guess_ka = 100;
ka = (low_guess_ka + high_guess_ka)/2;
if indx ~= 1
    tol = .1;
    T_max = t_max/2;
    while T_max < t_max - tol * t_max || T_max > t_max + tol * t_max
        if ka == 0 || ka == 100
            error('Ka solver hit the bounds. Go to extended1C_noMet.m to change low and high guess.')
        end
        ka = (low_guess_ka + high_guess_ka)/2;
        ode = @(t, y) [-ka*y(1); ...
            ka*F*y(1)/Vd - ke*y(2);...
            ke*y(2)*F_met*(Vd/Vd_met) - ke_met*y(3)];
        [T, Y] = ode45(ode, [0 1000], [Dose(1), 0, 0]);
        [argmax, maxVal] = max(Y(:,2));
        T_max = T(maxVal);
        if T_max < t_max - tol * t_max
            high_guess_ka = ka;
        elseif T_max > t_max + tol * t_max
            low_guess_ka = ka;
        end
    end
end

%IV infusion or patch
if indx == 1
    t = [];
    C = [[]];
    for i=1:num_doses-1
        sys = @(t, x) [-ke*x(1);...
               ke*F_met*(Vd/Vd_met)*x(1) - ke_met*x(2)];
        tspan = [dosing_times(i) dosing_times(i+1)];
        if i == 1
            y = [Dose(i)/Vd 0];
        else
            y = [C(end,1)+Dose(i)/Vd C(end,2)];
        end
        [temp_t, temp_C] = ode45(sys, tspan, y);
        t = [t; temp_t];
        C = [C; temp_C];
    end
    sys = @(t, x) [-ke*x(1);...
          ke*F_met*(Vd/Vd_met)*x(1) - ke_met*x(2)];
    tspan = [dosing_times(num_doses) range];
    y = [C(end,1)+Dose(num_doses)/Vd C(end,2)];
    [temp_t, temp_C] = ode45(sys, tspan, y);
    t = [t; temp_t];
    C = [C; temp_C];
else %Oral extended tablet
    t = [];
    C = [[]];
    for i=1:num_doses-1
        sys = @(t, x) [-ka*x(1); ...
            ka*F*x(1)/Vd - ke*x(2);...
            ke*F_met*(Vd/Vd_met)*x(2) - ke_met*x(3)];
        tspan = [dosing_times(i) dosing_times(i+1)];
        if i == 1
            y = [Dose(i) 0 0];
        else
            y = [(C(end,1)+Dose(i)) C(end,2) C(end,3)];
        end
        [temp_t, temp_C] = ode45(sys, tspan, y);
        t = [t; temp_t];
        C = [C; temp_C];
    end
    sys = @(t, x) [-ka*x(1); ...
            ka*F*x(1)/Vd - ke*x(2);...
            ke*F_met*(Vd/Vd_met)*x(2) - ke_met*x(3)];
    tspan = [dosing_times(num_doses) range];
    y = [C(end,1)+Dose(num_doses) C(end,2) C(end,3)];
    [temp_t, temp_C] = ode45(sys, tspan, y);
    t = [t; temp_t];
    C = [C; temp_C];
end

figure(1)
plot(t, C(:,end-1), 'b');
hold on; grid on;
plot(t, C(:,end), 'r');
xlabel('Time (hr)')
ylabel('Conc in central compartment (mg/L)');
switch indx
    case 1
        title('Drug Conc vs Time from IV Administration')
    case 2
        title('Drug Conc vs Time from Oral Administration')
    case 3
        title('Drug Conc vs Time from Intramusclar Administration')
    case 4
        title('Drug Conc vs Time from Subcutaneous Administration')
    case 5
        title('Drug Conc vs Time from Buccal Administration')
    case 6
        title('Drug Conc vs Time from Sublingual Administration')
    case 7
        title('Drug Conc vs Time from Lung Administration')
    case 8
        title('Drug Conc vs Time from Rectal Administration')
end
legend('Main Drug', 'Metabolite', 'Location', 'best')
xlim([0 range])
hold off;

figure(2)
plot(t, log(C(:,end-1)), 'b');
hold on; grid on;
plot(t, log(C(:,end)), 'r');
xlabel('Time (hr)')
ylabel('Log conc in central compartment (mg/L)');
switch indx
    case 1
        title('Drug Log Conc vs Time from IV Administration')
    case 2
        title('Drug Log Conc vs Time from Oral Administration')
    case 3
        title('Drug Log Conc vs Time from Intramusclar Administration')
    case 4
        title('Drug Log Conc vs Time from Subcutaneous Administration')
    case 5
        title('Drug Log Conc vs Time from Buccal Administration')
    case 6
        title('Drug Log Conc vs Time from Sublingual Administration')
    case 7
        title('Drug Log Conc vs Time from Lung Administration')
    case 8
        title('Drug Log Conc vs Time from Rectal Administration')
end
legend('Main Drug', 'Metabolite', 'Location', 'best')
xlim([0 range])
hold off;

end