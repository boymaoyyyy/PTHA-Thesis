from slip_generation.slip_sampler import StochasticSlipGenerator

sampler = StochasticSlipGenerator(magnitude=8.0,
    fault_params_csv='data/fault_geometry/heidarzadeh_2025_table2.csv')
print('slip_mean_m', sampler.slip_mean_m)
print('M0_target', sampler.M0_target)
slip = sampler.generate_rqmc_sample(0)
print('sample slip mean', slip.mean())
print('sample slip max', slip.max())
