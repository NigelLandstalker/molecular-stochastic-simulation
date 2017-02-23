#Reaction Simulation code for EE 5393 HW 1
#Author: Owen Hoffend

import math
import random
import operator as op
from functools import reduce

def k_choose_n(k, n):
	return math.factorial(k) / (math.factorial(n) * math.factorial(k - n))

def reaction_probs(reactant_descs, counts, k_values): #Returns the probability of a reaction happening
	#Reactant desc: tuple containing: ((r1 coefficient, r1 index), (r2 coefficient, r2 index), ..., (rn coefficient, rn index), k_value)
	#index gives the molecule type, while coefficient gives the number of molecules requied.
	alphas = []
	for index, desc in enumerate(reactant_descs):
		if can_fire(desc, counts):
			if type(desc[0]) == tuple:
				alphas.append(k_values[index] * reduce(op.mul, [k_choose_n(counts[desc[i][1]], desc[i][0]) for i in range(len(desc))]))
			else:
				alphas.append(k_values[index] * k_choose_n(counts[desc[1]], desc[0]))
		else:
			alphas.append(0)
	if sum(alphas) == 0:
		return 'end'
	return [alphas[i] / sum(alphas) for i in range(len(alphas))]

def can_fire(reactant_descs, counts): #Determines if a reaction can actually physically happen given the number of remaining molecules in the systen and the given reaction's chemical equation.
	can_fire_ = True
	if type(reactant_descs[0]) == tuple:
		for reactant in reactant_descs:
			if counts[reactant[1]] < reactant[0]:
				can_fire_ = False
	else:
		if counts[reactant_descs[1]] < reactant_descs[0]:
			can_fire_ = False
	return can_fire_

def stochastic_sim(init_counts, reaction_descs, iterations, end_conditions={}, return_intermediary_steps=False): #Runs a single stochastic simulation based on parameters
	#Reaction descs: ((reactant_descs), (product_descs), (k_values)). Tuple description of reaction. Reactant_descs is broken down into (coefficient, reactant_index).
	reactant_descs = [desc[0] for desc in reaction_descs]
	product_descs = [desc[1] for desc in reaction_descs]
	k_values = [desc[2] for desc in reaction_descs]
	counts = list(init_counts)
	if return_intermediary_steps:
		fired_probs = [] #Gives a list of the used
		all_counts = []

	for i in range(iterations):
		#print("Molecule counts for iteration %s: %s" % (iteration, counts))
		#probs = p1_probs(counts)
		probs = reaction_probs(reactant_descs, counts, k_values)
		#print("Reaction probs for iteration %s: %s" % (iteration, probs))
		if probs == 'end':
			print('Simulation finished after %s iterations' % i)
			break

		prob_sum = 0
		rand = random.uniform(0,1)
		for index, pb in enumerate(probs): #Generate a list of thresholded probabilities. For example: for probs 1/5, 2/5, 2/5, these would be: 1/5, 3/5, 5/5. Selection of a real number between 0 and 1 gives precisely one of the reactions
			prob_sum += pb
			if rand <= prob_sum:
				fired_reaction_index = index
				break

		fired_reactants = reactant_descs[fired_reaction_index]
		fired_products = product_descs[fired_reaction_index]
		if type(fired_reactants[0]) == tuple:
			for reactant in fired_reactants:
				counts[reactant[1]] -= reactant[0]
		else:
			counts[fired_reactants[1]] -= fired_reactants[0]

		if type(fired_products[0]) == tuple:
			for product in fired_products:
				counts[product[1]] += product[0]
		else:
			counts[fired_products[1]] += fired_products[0]

		if return_intermediary_steps:
			all_counts.append(list(counts))
			fired_probs.append(probs[fired_reaction_index])
	else:
		print("Simulation reached end of allowed iterations")

	if return_intermediary_steps:
		return (all_counts, fired_probs)
	else:
		return counts

def parse_reactions(reactions, counts={}, return_unique_molecules=False): #Just to make things easier
	#Reactions is a list of strings, each string corresponding to a reaction in the format: 02x1+03x2->01x2,k=1 where x is any alpha character
	#Reactions must completely desribe the chemical system desired! The parser assigns unique values to each molecular element.
	unique_molecules = {}
	parsed_equations = []
	for string in reactions:
		k_value = int(string.split(' ')[1][2:])
		string = string.split(' ')[0]

		reactants = string.split('->')[0].split('+')
		products = string.split('->')[1].split('+')

		def gen_tuple(molecule_list):
			list_of_tuples = []
			for m in molecule_list:
				if m[2:] not in unique_molecules:
					unique_molecules[m[2:]] = len(unique_molecules) #Assign new key value for a unique molecule
				list_of_tuples.append((int(m[:2]), unique_molecules[m[2:]]))
			return tuple(x for x in list_of_tuples)

		parsed_equations.append(((gen_tuple(reactants)), (gen_tuple(products)), k_value))

	if return_unique_molecules:
		return unique_molecules
	else:
		print("Parsed Equations: %s" % parsed_equations)
		if counts == {}:
			return parsed_equations
		else:
			parsed_counts = [0 for _ in range(len(unique_molecules))]
			for element in counts.keys():
				if element in unique_molecules:
					parsed_counts[unique_molecules[element]] = counts[element]
				else:
					print('Invalid molecule name in molecule counts')
			return [parsed_equations, parsed_counts]
			#Molecule count parsing happens here. Counts is a dict of keys with reactant names corresponding to the reactions above.

def p1_a_analyze_outcome(init_counts, reaction_descs, trial_count, iteration_count):
	#FIXME The conditions listed are TERMINATING CONDITIONS!
	all_outcomes = [stochastic_sim(init_counts, reaction_descs, iteration_count) for _ in range(trial_count)]

	C1_prob = len(list(filter(lambda x: x > 7, list(zip(*all_outcomes))[0]))) / trial_count
	C2_prob = len(list(filter(lambda x: x >= 8, list(zip(*all_outcomes))[1]))) / trial_count
	C3_prob = len(list(filter(lambda x: x < 3, list(zip(*all_outcomes))[2]))) / trial_count

	print("C1_prob: %s" % C1_prob)
	print("C2_prob: %s" % C2_prob)
	print("C3_prob: %s" % C3_prob)

def p1_b_analyze_outcome(init_counts, reaction_descs, trial_count, iteration_count):
	#Calculates a probability distribution for each variable that contains a state and its likelyhood of occurring
	prob_dist = {}
	encountered_paths = set()
	for iteration in range(trial_count):
		states = stochastic_sim(init_counts, reaction_descs, iteration_count, False, True) #States is structured: ([counts], [fired_probs])
		#All I really need to do is take the product of the correct index of probs and then make a set of all the possible ending configs?
		key = str(states[0][len(states[0]) - 1])
		value = reduce(op.mul, states[1])

		if key not in prob_dist:
			prob_dist[key] = value
		elif str(states[0]) not in encountered_paths:
			prob_dist[key] += value
		encountered_paths.add(str(states[0]))

	print(prob_dist)
	print("Number of unique results: %s" % len(prob_dist.values()))
	print("Searchspace Coverage: %s percent" % (sum(prob_dist.values()) * 100))
	return prob_dist
