import os
import pygad
import numpy as np
import random
import bisect
import argparse
from amino_to_dna_db import amino_acid_to_dna

import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run Pygad with specified parameters.")
    parser.add_argument("--num_generations", type=int, required=True, help="Number of generations for the GA")
    parser.add_argument("--sol_per_pop", type=int, required=True, help="Number of solutions per population")
    parser.add_argument("--num_parents_mating", type=int, required=True, help="Number of parents mating")
    parser.add_argument("--gene_space_high", type=str, required=True, help="Upper limit for gene space")
    parser.add_argument("--initial_pop_change_count", type=int, required=True, help="Initial population change count")
    parser.add_argument("--input_directory", type=str, required=True, help="Directory containing input sequences")
    parser.add_argument("--output_directory", type=str, required=True, help="Directory to save outputs")
    return parser.parse_args()


def is_dna_sequence(sequence):
    valid_dna_chars = {'A', 'T', 'C', 'G', 'N'}
    unique_chars = set(sequence)
    non_dna_chars = unique_chars - valid_dna_chars
    return len(non_dna_chars) <= 3           # 'A', 'T', 'C', 'G', 'N' 외에, 유니크한 문자열의 개수가 3개 이하여야 DNA 서열로 인식


def amino_acid_to_random_dna(amino_acid):
    if amino_acid in amino_acid_to_dna:
        return random.choice(amino_acid_to_dna[amino_acid])
    else:
        raise ValueError(f"Invalid amino acid: {amino_acid}")

def convert_amino_acid_sequence_to_dna(amino_acid_sequence):
    return ''.join([amino_acid_to_random_dna(aa) for aa in amino_acid_sequence])


def load_sequences_from_directory(directory_path):
    allowed_extensions = (".fa", ".txt")

    sequence_files = [f for f in os.listdir(directory_path) if f.endswith(allowed_extensions)]
    if not sequence_files:
        raise ValueError(f"No .fa or .txt files found in directory: {directory_path}")

    sequences = []
    sequence_ids = []
    for seq_file in sequence_files:
        with open(os.path.join(directory_path, seq_file), 'r') as f:
            lines = f.readlines()
            if not lines[0].startswith(">"):
                raise ValueError(f"Invalid FASTA format in {seq_file}")
            sequence_id = lines[0].strip()[1:]  # '>' 제거
            sequence_ids.append(sequence_id)
            sequence = ''.join(lines[1:]).replace('\n', '').strip()

            if is_dna_sequence(sequence):
                sequences.append(sequence)
            else:
                sequences.append(convert_amino_acid_sequence_to_dna(sequence))

    return sequences, sequence_ids

def map_char_to_number(char, custom_mapping, current_number):
    if char == 'A':
        return 1
    elif char == 'C':
        return 2
    elif char == 'G':
        return 3
    elif char == 'T':
        return 4
    else:
        if char not in custom_mapping:
            custom_mapping[char] = current_number[0]
            current_number[0] += 1
        return custom_mapping[char]

def fitness_function(solution, solution_idx, pad_sequence2, ga_instance):   # 각 숫자 비교해 같은 값을 가지면 +1, 아니면 0
    if not isinstance(pad_sequence2, (list, np.ndarray)):
        raise TypeError("pad_sequence2 must be a list or np.ndarray")
    
    fitness_score = sum([1 if solution[i] == pad_sequence2[i] else 0 for i in range(len(pad_sequence2))])

    if ga_instance.generations_completed == 0:  
        ga_instance.solutions = []  # 초기화
    ga_instance.solutions.append(fitness_score)  # solutions 배열에 적합도 값 추가
    

    return fitness_score

def make_initial_population(seq, sol_per_pop, custom_mapping, current_number, change_count):
    initial_population_ = []
    target_length = len(seq)

    change_count = min(change_count, target_length)

    for _ in range(sol_per_pop):
        result_str = list(seq)
        indices_to_change = random.sample(range(target_length), change_count)
        for idx in indices_to_change:
            current_char = result_str[idx]
            possible_chars = ["A", "T", "C", "G"]

            if current_char in possible_chars:
                possible_chars.remove(current_char)
                result_str[idx] = random.choice(possible_chars)
            else:
                print(f"Warning: {current_char} not in possible_chars. Skipping mutation.")

        numeric = [map_char_to_number(char, custom_mapping, current_number) for char in result_str]
        initial_population_.append(numeric)

    initial_population = np.array(initial_population_).reshape(sol_per_pop, -1)
    return initial_population



def weighted_sampler(seq, weights):
    totals = []
    for w in weights:
        if totals:
            totals.append(w + totals[-1])
        else:
            totals.append(w)
    return lambda: seq[bisect.bisect(totals, random.uniform(0, totals[-1]))]

def crossover_func(parents, offspring_size, ga_instance):
    offspring = []
    fitness_scores = ga_instance.last_generation_fitness
    sampler = weighted_sampler(ga_instance.population, fitness_scores)
    while len(offspring) != offspring_size[0]:
        parents = np.array([sampler() for _ in range(2)])
        parent1 = parents[0].copy()
        parent2 = parents[1].copy()
        random_split_point = np.random.choice(range(offspring_size[1]))
        result = np.concatenate((parent1[:random_split_point], parent2[random_split_point:]))
        offspring.append(result)
    return np.array(offspring)

def parent_selection_func(fitness, num_parents, ga_instance):
    num_population = ga_instance.population.shape[0]
    # select all population to parents
    selected_indices = np.arange(num_population)
    parents = ga_instance.population[selected_indices].copy()
    
    return parents, selected_indices

def mutation_func(offspring, ga_instance, pad_sequence2):
    base_characters = [1, 2, 3, 4, 5]
    
    new_offspring_fitness = []  

    for chromosome_idx in range(offspring.shape[0]):
        
        # Mutation: 2개의 랜덤한 값 변경 (기존 값과 다르게)
        random_gene_indices = np.random.choice(range(offspring.shape[1]), size=2, replace=False)
        
        for random_gene_idx in random_gene_indices:
            current_gene = offspring[chromosome_idx, random_gene_idx]
            available_characters = [c for c in base_characters if c != current_gene]  # 현재 값 제외
            new_gene = np.random.choice(available_characters)  
            offspring[chromosome_idx, random_gene_idx] = new_gene

        ###### Keep `5` at the end ######
        num_fives = np.sum(offspring[chromosome_idx] == 5)
        non_five_genes = offspring[chromosome_idx][offspring[chromosome_idx] != 5]
        offspring[chromosome_idx] = np.concatenate((non_five_genes, np.full(num_fives, 5)))

        ##################################### fitness #####################################
        offspring_fitness = fitness_function(offspring[chromosome_idx], chromosome_idx, pad_sequence2, ga_instance) 
        new_offspring_fitness.append(offspring_fitness)

    last_population = ga_instance.population
    last_population_fitness = ga_instance.last_generation_fitness.tolist()

    combined_list = [(value, index, list_num) 
                     for list_num, lst in enumerate([last_population_fitness, new_offspring_fitness]) 
                     for index, value in enumerate(lst)]
    combined_list.sort(reverse=True)
    
    top_elements = combined_list[:ga_instance.sol_per_pop - 1]
    
    result = []
    for value, index, list_num in top_elements:
        if list_num == 0:
            result.append(last_population[index].tolist())
        elif list_num == 1:
            result.append(offspring[index].tolist())
    
    final_offspring = np.array(result)
    return final_offspring


def on_generation(ga_instance):
    # 첫 번째 조건: 최대 세대 수에 도달했을 때
    if ga_instance.generations_completed >= ga_instance.num_generations:
        print(f"Stopping evolution at generation {ga_instance.generations_completed} (max generations reached).")
        ga_instance.generations_completed = ga_instance.num_generations  # 강제 종료

    # 두 번째 조건: stop_criteria에 맞는 조건이 만족됐을 때 (목표 길이를 찾았을 때)
    if ga_instance.solutions and 'reach_' + str(len(ga_instance.solutions)) in ga_instance.stop_criteria:
        print(f"Stopping evolution at generation {ga_instance.generations_completed} (target sequence reached).")
        ga_instance.generations_completed = ga_instance.num_generations  # 강제 종료

    # 로그 출력
    else:
        ga_instance.logger.info("Generation    = {generation}".format(generation=ga_instance.generations_completed))
        ga_instance.logger.info("Fitness_score = {fitness}".format(fitness=ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1]))

def run_ga(sequence, num_generations, sol_per_pop, output_path, sequence_id, num_parents_mating, gene_space_high, initial_pop_change_count):
    custom_mapping = {}
    current_number = [5]

    pad_sequence2 = [map_char_to_number(char, custom_mapping, current_number) for char in sequence]

    initial_population = make_initial_population(sequence, sol_per_pop, custom_mapping, current_number, initial_pop_change_count)

    ga_instance = pygad.GA(
        num_generations=num_generations,
        sol_per_pop=sol_per_pop,
        num_parents_mating=num_parents_mating,
        fitness_func=lambda ga, sol, sol_idx: fitness_function(sol, sol_idx, pad_sequence2, ga),
        num_genes=len(pad_sequence2),
        gene_space={"low": 1, "high": gene_space_high},
        crossover_type=crossover_func,
        mutation_type=lambda offspring, ga_instance: mutation_func(offspring, ga_instance, pad_sequence2),
        initial_population=initial_population,
        parent_selection_type=parent_selection_func,
        keep_parents=0,
        mutation_probability=1,
        gene_type=int,
        keep_elitism=1,
        stop_criteria=[f"reach_{len(sequence)}"],
        on_generation=on_generation
    )

    ga_instance.run()

    best_solution, best_solution_fitness, _ = ga_instance.best_solution()
    if best_solution_fitness == len(pad_sequence2):
        return True
    else:
        print(f"Target sequence {sequence_id} could not be reached through evolution.")
        with open(output_path, "a") as f:
            f.write(f"{sequence_id}\n")
        return False

def compare_sequences(directory_path, output_directory, num_generations, sol_per_pop, num_parents_mating, gene_space_high, initial_pop_change_count):
    sequences, sequence_ids = load_sequences_from_directory(directory_path)
    os.makedirs(output_directory, exist_ok=True)
    output_path = os.path.join(output_directory, "pygad_output.txt")

    all_reachable = True
    for idx, (sequence, sequence_id) in enumerate(zip(sequences, sequence_ids), start=1):
        remaining = len(sequences) - idx
        print(f"Processing sequence {sequence_id} ({remaining} remaining)...")
        
        # Use the passed `gene_space_high` parameter correctly
        current_gene_space_high = len(sequence) + 1 if gene_space_high is None else gene_space_high

        reachable = run_ga(
            sequence,
            num_generations,
            sol_per_pop,
            output_path,
            sequence_id,
            num_parents_mating,
            current_gene_space_high,
            initial_pop_change_count
        )
        if not reachable:
            all_reachable = False

    if all_reachable:
        print("All sequences could evolve to match themselves. No output file generated.")
    else:
        print(f"Results have been saved to {output_path}.")

def main():
    args = parse_arguments()

    NUM_GENERATIONS = args.num_generations
    SOL_PER_POP = args.sol_per_pop
    NUM_PARENTS_MATING = args.num_parents_mating
    GENE_SPACE_HIGH = args.gene_space_high
    INITIAL_POP_CHANGE_COUNT = args.initial_pop_change_count
    INPUT_DIRECTORY = args.input_directory
    OUTPUT_DIRECTORY = args.output_directory

    try:
        compare_sequences(
            INPUT_DIRECTORY,
            OUTPUT_DIRECTORY,
            NUM_GENERATIONS,
            SOL_PER_POP,
            NUM_PARENTS_MATING,
            GENE_SPACE_HIGH,
            INITIAL_POP_CHANGE_COUNT
        )
    except Exception as e:
        print(f"Error during execution: {e}")
main()


