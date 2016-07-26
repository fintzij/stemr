# Batch file to execute test_simulate_gillespie

library(batch)

for (sim_num in 1:2) {
        rbatch("./test_simulate_gillespie.R", seed = sim_num, sim_num = sim_num)
}

