import numpy
import msprime
import demesdraw
from matplotlib import pyplot as plt

# add demography
demography = msprime.Demography()
demography.add_population(name="N1", initial_size=100_000)
demography.add_population(name="N2", initial_size=100_000/3)

# instantaneous reduction of size 
#demography.add_population_parameters_change(population="EUR", time=100, initial_size=10_000)

demography.add_population(name="ANC", initial_size=7_000_000)
demography.add_population_split(time=5000, derived=["N1", "N2"], ancestral="ANC") # split 2k generations ago

# instatanous growth
#demography.add_population_parameters_change(time=4000, population="ANC", initial_size=100_000)

# instantaneous bottleneck
print(demography)

# Plot a schematic of the model
plt.clf()
demesdraw.tubes(demography.to_demes(), ax=plt.gca(), seed=1, log_time=True)
plt.show()

ts = msprime.sim_ancestry(
        {"N1": 10, "N2": 10}, 
        demography=demography, 
        recombination_rate=8.4e-9, # as in humans
        sequence_length=1_000,
        random_seed=1234)
print(ts)

# we can add mutations
mts = msprime.sim_mutations(ts, rate=3.5e-9, random_seed=1234)
print(mts.tables.sites)

# show SNPs
for variant in mts.variants():
    print(variant)

# visualise the haplotypes
samples = mts.samples()
for sample_id, h in zip(samples, mts.haplotypes(samples=samples)):
    pop = ts.node(sample_id).population
    print(f"Sample {sample_id:<2} ({ts.population(pop).metadata['name']:^5}): {h}")


# visualise the site frequency spectrum
plt.clf()
afs = mts.allele_frequency_spectrum()
plt.bar(range(mts.num_samples + 1), afs)
plt.title("Allele frequency spectrum")
plt.show()


# Define the samples between which Fst will be calculated
pop_id = {p.metadata["name"]: p.id for p in mts.populations()}
sample_sets=[mts.samples(pop_id["AFR"]), mts.samples(pop_id["EUR"])]

print(mts.Fst(sample_sets))


# Try more examples from https://tskit.dev/tutorials/popgen.html

######### IN LOOP

#### Now imposing distributions on each parameter

# T_split between 1000 and 8000 gens
tsplit_dist=np.random.uniform(low=1_000, high=8_000, size=8_000-1_000)
 
# N1
n1_size_dist=np.random.uniform(low=50_000, high=200_000, size=200_000-50_000)
plt.hist(n1_size_dist, bins=50, color='skyblue', alpha=0.7, rwidth=0.85)
plt.show()

# N2
n2_size_dist=n1_size_dist/30 ## Divide by 30
plt.hist(n2_size_dist, bins=50, color='skyblue', alpha=0.7, rwidth=0.85)
plt.show()

n2_size_dist=np.random.uniform(low=1666, high=6666, size=6666-1666) ## unsure about this
plt.hist(n2_size_dist, bins=50, color='skyblue', alpha=0.7, rwidth=0.85)
plt.show()

def repeat_simulations(mut, sample_sizes, length, reco, pop_size, num_simulations, seed=None):
    results = []
    for i in tqdm(range(num_simulations), desc="Running simulations"): 
        if seed is not None:
            np.random.seed(seed + i) 
            
            
        # Define demography
        # add demography
        demography = msprime.Demography()
        demography.add_population(name="N1", initial_size=100_000)
        demography.add_population(name="N2", initial_size=100_000/3)
        demography.add_population(name="ANC", initial_size=7_000_000)
        demography.add_population_split(time=5000, derived=["N1", "N2"], ancestral="ANC") # split 2k generations ago
        
        
        # Simulate 10 diploid samples under the coalescent with recombination on a 10kb region.
        ts = msprime.sim_ancestry(
            demography=demography,
            samples={"N1": sample_sizes, "N2": sample_sizes},
            recombination_rate=reco,
            sequence_length=length,
           # population_size=pop_size,
            random_seed=np.random.randint(99999999))
        
        # we can add mutations
        mutated_ts = msprime.sim_mutations(ts, rate=mut, random_seed=np.random.randint(99999999))

        diversity = mutated_ts.diversity()
        tajimas_d = mutated_ts.Tajimas_D()
        allele_frequency_spectrum = mutated_ts.allele_frequency_spectrum(polarised=True)
        results.append((mutated_ts, None, diversity, tajimas_d, allele_frequency_spectrum))
    return results



mut= 3.5e-9
length = 1_000
seed = 4710
reco = 8.4e-4
sample_sizes=10
pop_size = 1_000
num_simulations = 1

results = repeat_simulations(mut, sample_sizes, length, reco, pop_size,num_simulations, seed=seed)


