import re, yaml
import pandas as pd
import matplotlib.pyplot as plt
from configuration import *

# Generate figures using the results

print("Print Memory Footprint information")


def get_memory_footprint_info(filename):
    # Read the file
    with open(filename, "r") as file:
        content = file.read()
        # Use regular expressions to find the values
        dataset_memory_used_match = re.search(r"dataset-memory-used:\s*(\d+)", content)
        max_resident_set_size_match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", content)
        # Extract the values if found
        if dataset_memory_used_match and max_resident_set_size_match:
            dataset_memory_used = int(dataset_memory_used_match.group(1))
            max_resident_set_size = int(max_resident_set_size_match.group(1))
            algorithm_memory_used = (max_resident_set_size - dataset_memory_used) / (1024 * 1024)  # Convert to GB
            return algorithm_memory_used
        else:
            print("Error while reading the file")
            return 0


memory_footprint_data = {}
for dataset in datasets:
    memory_footprint_data[get_dataset_name(dataset)] = {}

for dataset in datasets:
    for algo in vtk_algorithms:
        algorithm_name = algorithms_names[algo]
        # TODO change to this
        # output_file = f"{memory_footprint_dir}/{get_dataset_name(dataset)}_{algorithm_name}.txt"
        output_file = f"{memory_footprint_dir}/{get_dataset_name(dataset)}_{algo}.txt"
        memory_footprint_data[get_dataset_name(dataset)][algorithm_name] = get_memory_footprint_info(output_file)

    # TODO Uncomment this
    # for algo in vtkm_algorithms:
    #     for hash_function in hash_functions[1:]:
    #         algorithm_name = algorithms_names[algo] + "-" + hash_function_names[hash_function]
    #         output_file = f"{memory_footprint_dir}/{get_dataset_name(dataset)}_{algorithm_name}.txt"
    #         memory_footprint_data[get_dataset_name(dataset)][algorithm_name] = get_memory_footprint_info(output_file)

# Convert the data to a pandas DataFrame
df_memory_fooprint = pd.DataFrame(memory_footprint_data)
df_memory_fooprint = df_memory_fooprint.transpose()
df_memory_fooprint.index.name = 'Dataset'
df_memory_fooprint.to_csv(f"{memory_footprint_dir}/memory-footprint.csv", index=True, header=True)
print(df_memory_fooprint)


# # Plot a line chart
# plt.figure(figsize=(12, 6))
# for algorithm in df_memory_fooprint.columns:
#     plt.plot(df_memory_fooprint.index, df_memory_fooprint[algorithm], marker='o', label=algorithm)
#     # Place algorithm names at the end of each line
#     # plt.text(df_memory_fooprint.index[-1], df_memory_fooprint[algorithm].values[-1], algorithm,
#     #          fontsize=10, va='center', ha='left')
#
# # Add labels and title
# # plt.title('Memory Footprint by Algorithm and Dataset', fontsize=14)
# plt.xlabel('Datasets', fontsize=12)
# plt.ylabel('Memory Footprint (gbytes)', fontsize=12)
# plt.xticks(rotation=45, ha='right')
# plt.legend(title='Algorithms', bbox_to_anchor=(1, 1), loc='upper left', fontsize=10)
# plt.grid()
# # Show the chart
# plt.tight_layout()
# plt.show()


def get_experiment_average_time_of_algorithm(experiment):
    algorithm_name = experiment['algorithm-name']
    hash_function = experiment['hash-name']
    algorithm_full_name = f"{algorithm_name}-{hash_function}" if hash_function != "None" else algorithm_name
    trials = experiment.get('trials', [])
    total_time = sum(trial['seconds-total'] for trial in trials)
    avg_time = total_time / len(trials) if trials else 0
    return algorithm_full_name, avg_time


def get_run_time_info(filename):
    # Read the yaml file
    with open(filename, "r") as file:
        content = yaml.load(file, Loader=yaml.FullLoader)[0]
        experiments = content["experiments"]
        # print(content)
        cpu_time_data_per_algorithm = {}
        for experiment in experiments:
            algorithm_full_name, avg_time = get_experiment_average_time_of_algorithm(experiment)
            cpu_time_data_per_algorithm[algorithm_full_name] = avg_time
            # print(f"Algorithm: {algorithm_full_name}, Average seconds-total: {avg_time:.4f}")
        return cpu_time_data_per_algorithm


print("Print CPU Time information")

cpu_time_data_1_threads_normal = {}
cpu_time_data_max_threads_normal = {}
cpu_time_data_max_threads_random = {}
for dataset in datasets:
    cpu_time_data_1_threads_normal[get_dataset_name(dataset)] = {}
    cpu_time_data_max_threads_normal[get_dataset_name(dataset)] = {}
    cpu_time_data_max_threads_random[get_dataset_name(dataset)] = {}

for dataset in datasets:
    cpu_time_data_1_threads_normal[get_dataset_name(dataset)] = get_run_time_info(
        f"{cpu_time_dir}/{get_dataset_name(dataset)}_1_threads_normal.yaml")
    cpu_time_data_max_threads_normal[get_dataset_name(dataset)] = get_run_time_info(
        f"{cpu_time_dir}/{get_dataset_name(dataset)}_{max_number_of_threads}_threads_normal.yaml")
    cpu_time_data_max_threads_random[get_dataset_name(dataset)] = get_run_time_info(
        f"{cpu_time_dir}/{get_dataset_name(dataset)}_{max_number_of_threads}_threads_random.yaml")

# Convert the data to a pandas DataFrame
df_1_threads_normal = pd.DataFrame(cpu_time_data_1_threads_normal)
df_1_threads_normal = df_1_threads_normal.transpose()
df_1_threads_normal.index.name = 'Dataset'
df_1_threads_normal.to_csv(f"{cpu_time_dir}/1_threads_normal.csv", index=True, header=True)

df_max_threads_normal = pd.DataFrame(cpu_time_data_max_threads_normal)
df_max_threads_normal = df_max_threads_normal.transpose()
df_max_threads_normal.index.name = 'Dataset'
df_max_threads_normal.to_csv(f"{cpu_time_dir}/{max_number_of_threads}_threads_normal.csv", index=True, header=True)

df_max_threads_random = pd.DataFrame(cpu_time_data_max_threads_random)
df_max_threads_random = df_max_threads_random.transpose()
df_max_threads_random.index.name = 'Dataset'
df_max_threads_random.to_csv(f"{cpu_time_dir}/{max_number_of_threads}_threads_random.csv", index=True, header=True)

cpu_time_dataframes = [df_1_threads_normal, df_max_threads_normal, df_max_threads_random]

for df in cpu_time_dataframes:
    print(df)
    # # Plot a line chart
    # plt.figure(figsize=(12, 6))
    # for algorithm in df.columns:
    #     plt.plot(df.index, df[algorithm], marker='o', label=algorithm)
    #     # Place algorithm names at the end of each line
    #     # plt.text(df.index[-1], df[algorithm].values[-1], algorithm,
    #     #          fontsize=10, va='center', ha='left')
    #
    # # Add labels and title
    # # plt.title('CPU time by Algorithm and Dataset', fontsize=14)
    # plt.xlabel('Datasets', fontsize=12)
    # plt.ylabel('Time (seconds)', fontsize=12)
    # plt.xticks(rotation=45, ha='right')
    # plt.legend(title='Algorithms', bbox_to_anchor=(1, 1), loc='upper left', fontsize=10)
    # plt.grid()
    # # Show the chart
    # plt.tight_layout()
    # plt.show()

print("Print hash distribution information")


def get_face_hash_distribution(filename):
    # Read the yaml file
    with open(filename, "r") as file:
        content = yaml.load(file, Loader=yaml.FullLoader)[0]
        experiments = content["experiments"]
        face_hash_distribution = experiments["face-hash-distribution"]
        return face_hash_distribution


for dataset in biggest_datasets:  # TODO change to datasets
    hash_distribution_file = f"{hash_performance_dir}/{get_dataset_name(dataset)}_hash_distribution.yaml"
    with open(hash_distribution_file, "r") as file:
        content = yaml.load(file, Loader=yaml.FullLoader)[0]
        experiments = content["experiments"]
        face_hash_distribution = experiments["face-hash-distribution"]
        print(f"Dataset: {get_dataset_name(dataset)}, Face Hash Distribution: {face_hash_distribution}")

print("Print cache misses information")


def get_cache_misses_info(filename):
    # Read the file content
    with open(filename, 'r') as file:
        text = file.read()
        # Use a regular expression to find the cache-misses value
        match = re.search(r'([\d,]+)\s+cache-misses', text)
        if match:
            return int(match.group(1).replace(",", ""))
        else:
            return 0


cache_misses_data = {}
for dataset in biggest_datasets:  # TODO change to datasets
    cache_misses_data[get_dataset_name(dataset)] = {}

for dataset in biggest_datasets:  # TODO change to datasets
    dataset_name = get_dataset_name(dataset)
    # cache misses for reading the dataset
    dataset_cache_misses_file = f"{hash_performance_dir}/{dataset_name}_cache_misses.txt"
    dataset_cache_misses = get_cache_misses_info(dataset_cache_misses_file)
    for algo in vtkm_algorithms:
        for hash_function in hash_functions[1:]:
            # TODO change this
            # algorithm_name = f"{algorithms_names[algo]}-{hash_function_names[hash_function]}"
            algorithm_name = f"{algo}_{str(hash_function)}"
            algorithm_cache_misses_file = f"{hash_performance_dir}/{dataset_name}_{algorithm_name}_cache_misses.txt"
            total_cache_misses = get_cache_misses_info(algorithm_cache_misses_file)
            algo_cache_misses = total_cache_misses - dataset_cache_misses
            cache_misses_data[dataset_name][algorithm_name] = algo_cache_misses
            # print(
            #     f"Dataset: {dataset_name}, Algorithm: {algorithm_name}, Cache Misses 1) DataSet: {dataset_cache_misses}, 2) Algo: {algo_cache_misses}")

# Convert the data to a pandas DataFrame
df_cache_misses = pd.DataFrame(cache_misses_data)
df_cache_misses = df_cache_misses.transpose()
df_cache_misses.index.name = 'Dataset'
df_cache_misses.to_csv(f"{hash_performance_dir}/cache-misses.csv", index=True, header=True)
print(df_cache_misses)

# # Plot a line chart
# plt.figure(figsize=(12, 6))
# for algorithm in df_cache_misses.columns:
#     plt.plot(df_cache_misses.index, df_cache_misses[algorithm], marker='o', label=algorithm)
#     # Place algorithm names at the end of each line
#     # plt.text(df_cache_misses.index[-1], df_cache_misses[algorithm].values[-1], algorithm,
#     #          fontsize=10, va='center', ha='left')
#
# # Add labels and title
# # plt.title('Cache Misses by Algorithm and Dataset', fontsize=14)
# plt.xlabel('Datasets', fontsize=12)
# plt.ylabel('Cache Misses', fontsize=12)
# plt.xticks(rotation=45, ha='right')
# plt.legend(title='Algorithms', bbox_to_anchor=(1, 1), loc='upper left', fontsize=10)
# plt.grid()
# # Show the chart
# plt.tight_layout()
# plt.show()

print("Print Parallel Efficiency information")
# TODO

print("Print GPU Time information")

# gpu_time_data_normal = {}
# gpu_time_data_random = {}
# for dataset in datasets:
#     gpu_time_data_normal[get_dataset_name(dataset)] = {}
#     gpu_time_data_random[get_dataset_name(dataset)] = {}
#
# for dataset in datasets:
#     gpu_time_data_normal[get_dataset_name(dataset)] = get_run_time_info(
#         f"{gpu_time_dir}/{get_dataset_name(dataset)}_normal.yaml")
#     gpu_time_data_random[get_dataset_name(dataset)] = get_run_time_info(
#         f"{gpu_time_dir}/{get_dataset_name(dataset)}_random.yaml")
#
# # Convert the data to a pandas DataFrame
# df_gpu_time_normal = pd.DataFrame(gpu_time_data_normal)
# df_gpu_time_normal = df_gpu_time_normal.transpose()
# df_gpu_time_normal.index.name = 'Dataset'
# df_gpu_time_normal.to_csv(f"{gpu_time_dir}/gpu_time_normal.csv", index=True, header=True)
#
# df_gpu_time_random = pd.DataFrame(gpu_time_data_random)
# df_gpu_time_random = df_gpu_time_random.transpose()
# df_gpu_time_random.index.name = 'Dataset'
# df_gpu_time_random.to_csv(f"{gpu_time_dir}/gpu_time_random.csv", index=True, header=True)
#
# gpu_time_dataframes = [df_gpu_time_normal, df_gpu_time_random]
#
# for df in gpu_time_dataframes:
#     print(df)
#     # # Plot a line chart
#     # plt.figure(figsize=(12, 6))
#     # for algorithm in df.columns:
#     #     plt.plot(df.index, df[algorithm], marker='o', label=algorithm)
#     #     # Place algorithm names at the end of each line
#     #     # plt.text(df.index[-1], df[algorithm].values[-1], algorithm,
#     #     #          fontsize=10, va='center', ha='left')
#     #
#     # # Add labels and title
#     # # plt.title('GPU time by Algorithm and Dataset', fontsize=14)
#     # plt.xlabel('Datasets', fontsize=12)
#     # plt.ylabel('Time (seconds)', fontsize=12)
#     # plt.xticks(rotation=45, ha='right')
#     # plt.legend(title='Algorithms', bbox_to_anchor=(1, 1), loc='upper left', fontsize=10)
#     # plt.grid()
#     # # Show the chart
#     # plt.tight_layout()
#     # plt.show()