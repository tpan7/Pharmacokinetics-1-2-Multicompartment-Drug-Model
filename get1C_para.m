function [parameters] = get1C_para(metabolite)
    prompt1 = 'What is the volume of distribution of the drug (L/kg)?: ';
    parameters(1) = input(prompt1)*70;
    
    prompt2 = 'What is the half-life of the drug from IV experiment (hr)?: ';
    parameters(2) = input(prompt2);
    
    prompt3 = 'What is the bioavailability of the drug?: ';
    parameters(3) = input(prompt3);
    if parameters(3) > 1 || parameters(3) < 0
        error('Invalid bioavailability');
    end
    
    prompt4 = 'What is the tmax of the drug (hr)?: ';
    parameters(4) = input(prompt4);
    
    prompt5 = 'How long do you want to model concentration (hr)?: ';
    parameters(5) = input(prompt5);
    
    if metabolite
        prompt6 = 'What is the volume of distribution of the metabolite (L/kg)?: ';
        parameters(6) = input(prompt6)*70;

        prompt7 = 'What is the half-life of the metabolite from IV experiment (hr)?: ';
        parameters(7) = input(prompt7);

        prompt8 = 'What is the fraction of drug is converted to metabolite?: ';
        parameters(8) = input(prompt8);
        if parameters(8) > 1 || parameters(8) < 0
            error('Invalid fraction converted');
        end
    end
end