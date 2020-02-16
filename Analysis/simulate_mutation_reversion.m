%% mutant_frequency = simulate_mutation_reversion(number_of_generations,mutation_rate,reversion_rate)
% simluate the effects of reversion rate on the mutant frequency
% LBC Febuary 2020
function mutant_frequency = simulate_mutation_reversion(number_of_generations,mutation_rate,reversion_rate)
    WTfreq = 1 ;
    mutant_frequency = NaN( 1 , number_of_generations );
    mutant_frequency(1) = 1-WTfreq ;
    for I = 2:number_of_generations
        mutant_frequency(I) = mutant_frequency(I-1)  + ...
        (1-mutant_frequency(I-1)) * mutation_rate - reversion_rate*mutant_frequency(I-1) ;
    end

end