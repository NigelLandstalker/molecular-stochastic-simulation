#Reaction Simulation code for EE 5393 HW 1
#Author: Owen Hoffend

import math
import random
from functools import reduce

def k_choose_n(k, n):
		return math.factorial(k) / (math.factorial(n) * math.factorial(k - n))

def reaction_probs(reactant_descs, counts, k_values): #Returns the probability of a reaction happening
	#Reactant desc: tuple containing: ((r1 coefficient, r1 index), (r2 coefficient, r2 index), ..., (rn coefficient, rn index), k_value)
	#index gives the molecule type, while coefficient gives the number of molecules requied.
	alphas = []
	for index, desc in enumerate(reactant_descs):
		if can_fire(desc, counts):
			alphas.append(k_values[index] * reduce(lambda x, y: x * y, [k_choose_n(counts[desc[index][1]], desc[index][0]) for index in range(len(desc))]))
		else:
			alphas.append(0)
	if sum(alphas) == 0:
		return 'end'
	return [alphas[i] / sum(alphas) for i in range(len(alphas))]

def p1_probs(counts): #For testing purposes: using the probability equations given in the problem. Should be equivalent to reaction_probs!
	x1, x2, x3 = counts[0], counts[1], counts[2]
	denominator = (0.5 * x1 * (x1 - 1) * x2) + (x1 * x3 * (x3 - 1)) + (3 * x2 * x3)
	p1 = (0.5 * x1 * (x1 - 1) * x2) / denominator
	p2 = (x1 * x3 * (x3 - 1)) / denominator
	p3 = (3 * x2 * x3) / denominator
	return [p1, p2, p3]

def can_fire(reactant_descs, counts): #Determines if a reaction can actually physically happen given the number of remaining molecules in the systen and the given reaction's chemical equation.
	can_fire_ = True
	if type(reactant_descs) == tuple:
		for reactant in reactant_descs:
			if counts[reactant[1]] < reactant[0]:
				can_fire_ = False
	else:
		if counts[reactant_descs[1]] < reactant_descs[0]:
			can_fire_ = False
	return can_fire_

def stochastic_sim(init_counts, reaction_descs, iterations, halt_on_low_reactants=False, return_intermediary_steps=False): #Runs a single stochastic simulation based on parameters
	#Reaction descs: ((reactant_descs), (product_descs), (k_values)). Tuple description of reaction. Reactant_descs is broken down into (coefficient, reactant_index).
	reactant_descs = [desc[0] for desc in reaction_descs]
	product_descs = [desc[1] for desc in reaction_descs]
	k_values = [desc[2] for desc in reaction_descs]
	counts = list(init_counts)
	if return_intermediary_steps:
		fired_probs = [] #Gives a list of the used
		all_counts = []

	for _ in range(iterations):
		#if(0 in counts): break
		#print("Molecule counts for iteration %s: %s" % (iteration, counts))
		#probs = p1_probs(counts)
		probs = reaction_probs(reactant_descs, counts, k_values) #Probs should actually be a list containing the reaction index for the reaction as well as the probability of it ocurring
		#print("Reaction probs for iteration %s: %s" % (iteration, probs))
		if probs == 'end' or halt_on_low_reactants and False in [can_fire(desc, counts) for desc in reactant_descs]:
			break

		prob_sum = 0
		rand = random.uniform(0,1)
		for index, pb in enumerate(probs): #Generate a list of thresholded probabilities. For example: for probs 1/5, 2/5, 2/5, these would be: 1/5, 3/5, 5/5. Selection of a real number between 0 and 1 gives precisely one of the reactions
			prob_sum += pb
			if rand <= prob_sum:
				fired_reaction_index = index
				break

		#print("Reaction of index %s fired" % fired_reaction_index)
		fired_reactants = reactant_descs[fired_reaction_index]
		fired_products = product_descs[fired_reaction_index]
		for reactant in fired_reactants:
			counts[reactant[1]] -= reactant[0]

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

def p1_a_analyze_outcome(init_counts, reaction_descs, trial_count, iteration_count, halt_on_low_reactants):
	all_outcomes = [stochastic_sim(init_counts, reaction_descs, iteration_count, halt_on_low_reactants) for _ in range(trial_count)]

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
		value = reduce(lambda x, y: x * y, states[1])

		if key not in prob_dist:
			prob_dist[key] = value
		elif str(states[0]) not in encountered_paths:
			prob_dist[key] += value
		encountered_paths.add(str(states[0]))

	print(prob_dist)
	print("Number of unique results: %s" % len(prob_dist.values()))
	print("Searchspace Coverage: %s percent" % (sum(prob_dist.values()) * 100))
	return prob_dist

#Reactions for problem 1:
r1 = (((2,0), (1,1)),(4, 2), 1)
r2 = (((1,0), (2,2)),(3, 1), 2)
r3 = (((1,1), (1,2)),(2, 0), 3)
p1_descs = [r1, r2, r3]
p1_a_counts = [5, 5, 5]
#print(stochastic_sim(p1_counts, p1_descs, 1000))
#p1_a_analyze_outcome(p1_a_counts, p1_descs, 10000, 1000, True)

p1_b_counts = [6, 6, 6]
p1_b_analyze_outcome(p1_b_counts, p1_descs, 1000000, 5)

#Trial Results:
#Organic simulation:
#Command used: p1_a_analyze_outcome(p1_counts, p1_descs, number of trials, number of reactions per trial, True/False)
#1000 trials at 10 reactions per trial:
#C1: 0.16, C2: 0.32, C3: 0.32

#1000 trials at 100 reactions per trial:
#C1: 0.258, C2: 1.0, C3: 0.007

#1000 trials at 1000 reactions per trial:
#C1: 0.949, C2: 1.0, C3: 0.007

#Halting simulation at first instance of an exhausted reactant (first time any reaction can't happen):
#10000 Trials:
#C1: 0.191, C2: 0.78, C3: 0.541

#{'[9, 5, 5]': 0.09527239811349611, '[6, 9, 4]': 0.39955463479823056, '[13, 5, 0]': 0.0008974886503640765, '[1, 5, 15]': 0.004308255482897824, '[2, 9, 9]': 0.18325290571277889, '[0, 1, 21]': 7.807544039428096e-07, '[5, 5, 10]': 0.13499120346475102, '[16, 1, 1]': 3.873107321383183e-05, '[0, 17, 2]': 0.003281885587719334, '[4, 1, 16]': 0.0004894563339882729, '[12, 1, 6]': 0.0015374184639134947, '[3, 13, 3]': 0.17339259949008246, '[8, 1, 11]': 0.002981149841492497}
#Number of unique results: 13
#Searchspace Coverage: 99.99989077673322 percent
