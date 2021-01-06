function extended2C_noMet(parameters, route)
%Which route of administration
%1 for Intravenous, 2 for Oral, 3 for Intramusclar, 4 for Subcutaneous
%5 for Buccal, 6 for Sublingual, 7 for Lung, 8 for Rectal, 9 for
%transdermal
[temp, indx] = max(route);

%Blood flow rate
Q = 24;
%The Volume of distribution of first compartment
Vd1 = parameters(1);
%The Volume of distribution of 2nd compartment
Vd2 = parameters(2);
%The half life
t_half = parameters(3);
%The bioavailability factor
F = parameters(4);
if indx == 1
    F = 1;
end
%The time where max concentration occurs
t_max = parameters(5);
%The fraction unbound of drug
fu = parameters(6);
%The end of the plotting
range = parameters(7);

prompt_dose = 'What is the total dose of the drug (mg)?: ';
Dose = input(prompt_dose);

prompt_time = 'Over what length of time is the drug being administered (hr)? In other words, what is the extinction time of administration?: ';
t_extinct = input(prompt_time);

%Define IV experiment to solve for ke
ke_ode = @(ke, t, x) [(Q/Vd1)*(x(2)-fu*x(1)) - ke*fu*x(1);...
        (Q/Vd2)*(fu*x(1) - x(2))];
ke = ke_solve2C(parameters, ke_ode, [Dose/Vd1 0]);

%To solve for ka
if indx == 2
    ka_ode = @(ka, t, x) [-ka*x(1);...
            ka*F*x(1)/Vd1 + (Q/Vd1)*(x(3)-fu*x(2)) - ke*fu*x(2);...
            (Q/Vd2)*(fu*x(2) - x(3))];
    ka = ka_solve2C(parameters, ka_ode, [Dose 0 0]);
end


%IV infusion or patch
if indx == 1 || indx == 9
    %During infusion or patch
    sys1 = @(t, x) [(Dose/Vd1)*F/t_extinct + (Q/Vd1)*(x(2)-fu*x(1)) - ke*fu*x(1);...
        (Q/Vd2)*(fu*x(1) - x(2))];
    tspan1 = [0 t_extinct];
    y1 = [0 0];
    [t1, C1] = ode45(sys1, tspan1, y1);
    
    if range > t_extinct
        sys2 = @(t, x) [(Q/Vd1)*(x(2)-fu*x(1)) - ke*fu*x(1);...
        (Q/Vd2)*(fu*x(1) - x(2))];
        tspan2 = [t_extinct range];
        y2 = [C1(end,1) C1(end,2)];
        [t2, C2] = ode45(sys2, tspan2, y2);
        
        t1 = [t1; t2];
        C1 = [C1; C2];
    end
else %Oral extended tablet
    sys1 = @(t, x) [Dose/t_extinct-ka*x(1);...
            ka*F*x(1)/Vd1 + (Q/Vd1)*(x(3)-fu*x(2)) - ke*fu*x(2);...
            (Q/Vd2)*(fu*x(2) - x(3))];
    tspan1 = [0 t_extinct];
    y1 = [0 0 0];
    [t1, C1] = ode45(sys1, tspan1, y1);
    
    if range > t_extinct
        sys2 = @(t, x) [-ka*x(1);...
            ka*F*x(1)/Vd1 + (Q/Vd1)*(x(3)-fu*x(2)) - ke*fu*x(2);...
            (Q/Vd2)*(fu*x(2) - x(3))];
        tspan2 = [t_extinct range];
        y2 = [C1(end, 1) C1(end, 2) C1(end, 3)];
        [t2, C2] = ode45(sys2, tspan2, y2);
        
        t1 = [t1; t2];
        C1 = [C1; C2];
    end
end

figure(1)
plot(t1, C1(:,end-1), 'b');
hold on; grid on;
plot(t1, C1(:,end), 'r');
xlabel('Time (hr)')
ylabel('Conc of drug in central compartment (mg/L)');
switch indx
    case 1
        title('Drug Conc vs Time from IV Infusion')
    case 2
        title('Drug Conc vs Time from Oral Extended Tablet')
    case 9
        title('Drug Conc vs Time from Transdermal Patch')
end
xlim([0 range])
legend('Central Compartment', '2nd Compartment', 'Location', 'best')
hold off;

figure(2)
plot(t1, log(C1(:,end-1)), 'b')
hold on; grid on;
plot(t1, log(C1(:,end)), 'r')
xlabel('Time (hr)')
ylabel('Log conc of drug in central compartment (mg/L)');
switch indx
    case 1
        title('Drug Log Conc vs Time from IV Infusion')
    case 2
        title('Drug Log Conc vs Time from Oral Extended Tablet')
    case 9
        title('Drug Log Conc vs Time from Transdermal Patch')
end
xlim([0 range])
legend('Central Compartment', '2nd Compartment', 'Location', 'best')
hold off;

end