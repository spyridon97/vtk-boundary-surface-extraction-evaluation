from configuration import *

# Create directories

# Memory footprint information
print("Get memory footprint information")
os.makedirs(memory_footprint_dir, exist_ok=True)
for dataset in datasets:
    for algo in vtk_algorithms:
        run_command(f"{memory_evaluator} {executable} -i {dataset} -d TBB -t 1 {algo} -n 0",
                    f"{memory_footprint_dir}/{get_dataset_name(dataset)}_{algorithms_names[algo]}.txt")

    for algo in vtkm_algorithms:
        for hash_function in hash_functions[1:]:
            run_command(f"{memory_evaluator} {executable} -i {dataset} -d TBB -t 1 {algo} -f {hash_function} -n 0",
                        f"{memory_footprint_dir}/{get_dataset_name(dataset)}_{algorithms_names[algo]}-{hash_function_names[hash_function]}.txt")
# CPU time information
print("Get CPU time information")
os.makedirs(cpu_time_dir, exist_ok=True)
for dataset in datasets:
    run_command(f"{executable} -i {dataset} -d TBB -t 1 {algorithms_joined} -n {iterations}",
                f"{cpu_time_dir}/{get_dataset_name(dataset)}_1_threads_normal.yaml")
    run_command(
        f"{executable} -i {dataset} -d TBB -t {max_number_of_threads} {parallel_algorithms_joined} -n {iterations}",
        f"{cpu_time_dir}/{get_dataset_name(dataset)}_{max_number_of_threads}_threads_normal.yaml")
    run_command(
        f"{executable} -i {dataset} -d TBB -t {max_number_of_threads} {parallel_algorithms_joined} -n {iterations} -r",
        f"{cpu_time_dir}/{get_dataset_name(dataset)}_{max_number_of_threads}_threads_random.yaml")

# Hash distribution information
print("Get hash distribution information")
os.makedirs(hash_performance_dir, exist_ok=True)
for dataset in datasets:
    run_command(f"{executable} -i {dataset} -d TBB -t 1 --hash-distribution",
                f"{hash_performance_dir}/{get_dataset_name(dataset)}_hash_distribution.yaml")

# Cache misses information
print("Get cache misses information")
for dataset in datasets:
    # cache misses for reading the dataset
    run_command(f"{cache_evaluator} {executable} -i {dataset} -d TBB -n 0",
                f"{hash_performance_dir}/{get_dataset_name(dataset)}_cache_misses.txt")
    # cache misses for each algorithm with each hash function
    for algo in vtkm_algorithms:
        for hash_function in hash_functions[1:]:
            run_command(f"{cache_evaluator} {executable} -i {dataset} -d TBB -t 1 {algo} -f {hash_function} -n 0",
                        f"{hash_performance_dir}/{get_dataset_name(dataset)}_{algorithms_names[algo]}-{hash_function_names[hash_function]}_cache_misses.txt")

# Parallel efficiency information
print("Get time information for parallel algorithms")
os.makedirs(parallel_efficiency_dir, exist_ok=True)
for dataset in biggest_datasets:
    for power in range(0, max_number_of_threads_power_of_2 + 1):
        threads = int(math.pow(2, power))
        run_command(f"{executable} -i {dataset} -d TBB -t {threads} {parallel_algorithms_joined} -n {iterations}",
                    f"{parallel_efficiency_dir}/{get_dataset_name(dataset)}_{threads}_threads.yaml")

# GPU time information
print("Get GPU time information")
os.makedirs(gpu_time_dir, exist_ok=True)
for dataset in datasets:
    run_command(f"{executable} -i {dataset} -d KOKKOS {vtkm_algorithms_joined} -n {iterations}",
                f"{gpu_time_dir}/{get_dataset_name(dataset)}_normal.yaml")
    run_command(f"{executable} -i {dataset} -d KOKKOS {vtkm_algorithms_joined} -n {iterations} -r",
                f"{gpu_time_dir}/{get_dataset_name(dataset)}_random.yaml")

print(f"All evaluations completed. Check results in {results_dir}/.")
