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
			#r1_molecule_count = counts[desc[0][1]]
			#r2_molecule_count = counts[desc[1][1]]
			#alphas.append(k_values[index] * k_choose_n(r1_molecule_count, desc[0][0]) * k_choose_n(r2_molecule_count, desc[1][0]))
			alphas.append(k_values[index] * reduce(lambda x, y: x * y, [k_choose_n(counts[desc[index][1]], desc[index][0]) for index in range(len(desc))]))

	return [alphas[i] / sum(alphas) for i in range(len(alphas))]

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

def stochastic_sim(init_counts, reaction_descs, iterations):
	#Reaction descs: ((reactant_descs), (product_descs), (k_values))
	reactant_descs = [desc[0] for desc in reaction_descs]
	product_descs = [desc[1] for desc in reaction_descs]
	k_values = [desc[2] for desc in reaction_descs]
	counts = init_counts

	for iteration in range(iterations):
		probs = reaction_probs(reactant_descs, counts, k_values) #Probs should actually be a list containing the reaction index for the reaction as well as the probability of it ocurring

		print("Molecule counts for iteration %s: %s" % (iteration, counts))
		print("Reaction probs for iteration %s: %s" % (iteration, probs))

		prob_sum = 0
		rand = random.uniform(0,1)
		for index, pb in enumerate(probs): #Generate a list of thresholded probabilities. For example: for probs 1/5, 2/5, 2/5, these would be: 1/5, 3/5, 5/5. Selection of a real number between 0 and 1 gives precisely one of the reactions
			prob_sum += pb
			if rand <= prob_sum:
				fired_reaction_index = index
				break

		fired_reactants = reactant_descs[fired_reaction_index]
		if can_fire(fired_reactants, counts):
			fired_products = product_descs[fired_reaction_index]
			for reactant in fired_reactants:
				counts[reactant[1]] -= reactant[0]

			if type(fired_products[0]) == tuple:
				for product in fired_products:
					counts[product[1]] += product[0]
			else:
				counts[fired_products[1]] += fired_products[0]
	return counts

#Reactions for problem 1:
r1 = (((2,0), (1,1)),(4, 2), 1)
r2 = (((1,0), (2,2)),(3, 1), 2)
r3 = (((1,1), (1,2)),(2, 0), 3)
p1_descs = [r1, r2, r3]
p1_counts = [3, 3, 3]
print(stochastic_sim(p1_counts, p1_descs, 2))

#NEXT TASKS: Implement a "can_fire" function within the stochastic_sim function, since otherwise reactants will cause the products to go negative soon enough.
