from slip_generation.slip_sampler import StochasticSlipGenerator

g = StochasticSlipGenerator(8.0, 'data/fault_geometry/heidarzadeh_2025_table2.csv')
print('L_km', g.L_km, 'W_km', g.W_km, 'slip_mean', g.slip_mean_m)
print('n_along', g.n_along, 'n_down', g.n_down, 'N', g.N)
