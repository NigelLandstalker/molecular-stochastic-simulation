#Reaction Simulation code for EE 5393 HW 1
#Author: Owen Hoffend

import math
import random
import operator as op
from functools import reduce

def k_choose_n(k, n):
	#This function takes two arguments k and n and returns k choose n.
	return math.factorial(k) / (math.factorial(n) * math.factorial(k - n))

def reaction_probs(reactant_descs, counts, k_values):
	#Returns the probability of a reaction happening based on the reaction descriptions, the number of molecules of each reactant, and their respective k-values.
	#Reactant desc: tuple containing: ((r1 coefficient, r1 index), (r2 coefficient, r2 index), ..., (rn coefficient, rn index))
	#index gives the molecule type, while coefficient gives the number of molecules requied.
	#Function reaturns a list of probabilities between 0 and 1
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

def can_fire(reactant_descs, counts):
	#Determines if a reaction can actually physically happen given the number of remaining molecules in the systen and the given reaction's chemical equation.
	can_fire_ = True
	if type(reactant_descs[0]) == tuple:
		for reactant in reactant_descs:
			if counts[reactant[1]] < reactant[0]:
				can_fire_ = False
	else:
		if counts[reactant_descs[1]] < reactant_descs[0]:
			can_fire_ = False
	return can_fire_

def stochastic_sim(init_counts, reaction_descs, iterations, end_conditions=[], return_intermediary_steps=False):
	#Runs a single stochastic simulation based on parameters
	#Reaction descs: ((reactant_descs), (product_descs), (k_values)). Tuple description of reaction. Reactant_descs is broken down into (coefficient, reactant_index).
	reactant_descs = [desc[0] for desc in reaction_descs]
	product_descs = [desc[1] for desc in reaction_descs]
	k_values = [desc[2] for desc in reaction_descs]
	counts = list(init_counts)
	if return_intermediary_steps:
		fired_probs = [] #Gives a list of the used
		all_counts = []

	last_num = 0
	last_num_freq = 0
	print_outs = []
	for i in range(iterations):
		#print("Molecule counts for iteration %s: %s" % (iteration, counts))
		#probs = p1_probs(counts)
		probs = reaction_probs(reactant_descs, counts, k_values)
		#print("Reaction probs for iteration %s: %s" % (iteration, probs))
		ending_state = [function(counts) for function in end_conditions]
		if probs == 'end' or True in ending_state:
			return [ending_state, i]

		#Determine which reaction to fire:
		prob_sum = 0
		rand = random.uniform(0,1)
		for index, pb in enumerate(probs): #Generate a list of thresholded probabilities. For example: for probs 1/5, 2/5, 2/5, these would be: 1/5, 3/5, 5/5. Selection of a real number between 0 and 1 gives precisely one of the reactions
			prob_sum += pb
			if rand <= prob_sum:
				fired_reaction_index = index
				break

		fired_reactants = reactant_descs[fired_reaction_index]
		fired_products = product_descs[fired_reaction_index]

		#Print some stuff about the reactions (used for collatz conjecture)
		#if counts[8] > 0:
			#if counts[8] == last_num:
				#last_num_freq += 1
			#else:
				#last_num_freq = 0

			#if last_num_freq > 8:
				#if counts[8] not in print_outs:
					#print_outs.append(counts[8])
					#print(counts[8])
				#last_num_freq = 0
			#last_num = counts[8]

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

	if return_intermediary_steps:
		return (all_counts, fired_probs)
	else:
		return counts

def parse_reactions(reactions, counts={}, return_unique_molecules=False): #Just to make things easier
	#Reactions is a list of strings, each string corresponding to a reaction in the format: 02x1+03x2->01x2,k=1 where x is any alpha character
	#Reactions must completely desribe the chemical system desired! The parser assigns unique values to each molecular element.
	#The parser returns a list of parsed equations that can be used by the functions above.

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
		print("Successfully parsed equations")
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
