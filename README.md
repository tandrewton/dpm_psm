# dpm_psm
minimal example for 2024 PSM paper

# Notes
To use this directory, please follow the following instructions to run a minimal simulation.
Simulation processing code will not be made available in general, please refer to the output streams in the cpp files to format the data in your preferred manner.

# Compilation command:
g++ -O3 --std=c++11 -g -I src main/cell/psm2D.cpp src/dpm.cpp src/cell.cpp -o main/cell/psm2D.o

# Example of a script that will generate an array of run commands for various parameter choices

att_arr=(0.05)

att2_arr=(0.05)

t_stress_arr=(10000.0)

v0=0.1

phi_arr=(0.8)

tau_abp=1.0

gamma_arr=(0)

kon_arr=(1.0)

koff_arr=(0.5)

kl=0.5

ka=(5.0)

kb=0.01

calcMinPos=1

for att in ${att_arr[@]}; do
  for att2 in ${att2_arr[@]}; do
    for phi in ${phi_arr[@]}; do
      for t_stress in ${t_stress_arr[@]}; do
        for gamma in ${gamma_arr[@]}; do
          for k_on in ${kon_arr[@]}; do
            for k_off in ${koff_arr[@]}; do
              k_ecm=$att2
              echo "./main/cell/psm2D.o   12  16 1.0 $phi $kl $ka $kb $att $att2 $t_stress    $v0    $tau_abp  $gamma $k_on $k_off $k_ecm $calcMinPos 1    50 testa_"$att"_a2_"$att2"_tm_"$t_stress"_p_"$phi"_t_"$tau_abp"_gamma_"$gamma"_k_on_"$k_on"_k_off_"$k_off"_k_ecm_"$k_ecm
            done
          done
        done
      done
    done
  done
done

# an actual run command for a single simulation:


./main/cell/psm2D.o   12  16 1.0 0.8 0.5 5.0 0.01 0.05 0.05 10000.0    0.1    1.0  0 1.0 0.5 0.05 1 1   
 50 testa_0.05_a2_0.05_tm_10000.0_p_0.8_t_1.0_gamma_0_k_on_1.0_k_off_0.5_k_ecm_0.05
