#Molecular Logic Elements for EE 5393 HW 1
#Uses reaction_simulation.py
#Author: Owen Hoffend

from reaction_simulation import parse_reactions
from reaction_simulation import stochastic_sim

import operator as op
from functools import reduce
import ast

SLOW = '1'
MEDIUM = '100'
FAST = '10000'
VERY_FAST = '1000000'

#def call_reaction(reaction, counts, iterations):
#	sim_input = parse_reactions(reaction, counts)
#	output = stochastic_sim(sim_input[1], sim_input[0], iterations)
#	return output

def statistical_call_reaction(reaction, counts, iterations, trials):
	print("Initiating statistical reaction simulation of %s trials, with %s iterations per trial." % (trials, iterations))
	sim_input = parse_reactions(reaction, counts)
	outputs = [stochastic_sim(sim_input[1], sim_input[0], iterations) for _ in range(trials)]
	averages = [(sum(counts) / float(trials)) for counts in list(zip(*outputs))]
	print('Simulation finished')
	return associate_reactants(reaction, averages)

def associate_reactants(reaction, outputs): #Associates the reactants printed in the output back with their respective molecule string.
#This function only exists because of my bad code that I don't want to change.
	unique_molecules = parse_reactions(reaction, return_unique_molecules=True)
	new_output = [key + ": " + str(outputs[unique_molecules[key]]) for key in unique_molecules]
	return new_output

def abs_indicator(molecule):
	return [
		"00o->01" + molecule + "_ab' k=" + SLOW,
		"01" + molecule + "+01" + molecule + "_ab'->01" + molecule + " k=" + FAST,
		"02" + molecule + "_ab'->01" + molecule + "_ab'" + " k=" + FAST,
		#buffer reactions:
		"01" + molecule + "_ab'->01" + molecule+ "_ab k=" + SLOW,
		"01" + molecule + "+01" + molecule + "_ab->01" + molecule + " k=" + FAST,
		"02" + molecule + "_ab->01" + molecule + "_ab k=" + FAST
	]

def p1_a_analyze_outcome(trial_count, iteration_count):
	end_conditions = [
		lambda x: x[0] > 7,
		lambda x: x[1] >= 8,
		lambda x: x[2] < 3
	]
	sim_input = parse_reactions(part1_equations, {'x1':5,'x2':5,'x3':5})
	outputs = [stochastic_sim(sim_input[1], sim_input[0], iteration_count, end_conditions=end_conditions) for _ in range(trial_count)]

	conditions = [output[0] for output in outputs]
	halt_iterations = [output[1] for output in outputs]

	average_halt_iterations = sum(halt_iterations) / float(len(halt_iterations))
	condition_averages = [sum(map(lambda x: 1 if x else 0, cond)) / float(trial_count) for cond in list(zip(*conditions))]

	print("Condition probabilities: x1 > 7: %s, x2 >= 8: %s, x3 < 3: %s" % tuple(condition_averages))
	print("Average number of iterations before halting: %s" % average_halt_iterations)
	return (condition_averages, average_halt_iterations)

def p1_b_analyze_outcome(trial_count):
	sim_input = parse_reactions(part1_equations, {'x1':6,'x2':6,'x3':6})
	#Calculates a probability distribution for each variable that contains a state and its likelyhood of occurring
	prob_dist = {}
	encountered_paths = set()
	for iteration in range(trial_count):
		states = stochastic_sim(sim_input[1], sim_input[0], 5, return_intermediary_steps=True) #States is structured: ([counts], [fired_probs])
		key = str(states[0][len(states[0]) - 1])
		value = reduce(op.mul, states[1])

		if key not in prob_dist:
			prob_dist[key] = value
		elif str(states[0]) not in encountered_paths:
			prob_dist[key] += value
		encountered_paths.add(str(states[0]))

	end_states = [(ast.literal_eval(key), prob_dist[key]) for key in prob_dist]
	num_prob_dist = {}
	for i in range(3):
		elements = {}
		for element in end_states:
			if element[0][i] in elements:
				elements[element[0][i]] += element[1]
			else:
				elements[element[0][i]] = element[1]
		num_prob_dist['x' + str(i + 1)] = elements

		print('x' + str(i + 1) + " Probability distribution:")
		for key, value in sorted(elements.items()):
			print("{} : {}".format(key, value))

	print("Output states probability distribution: %s" % prob_dist)
	print("Number of unique results: %s" % len(prob_dist.values()))
	print("Searchspace Coverage: %s percent" % (sum(prob_dist.values()) * 100))

	#individual number probability distribution:

	return (prob_dist, num_prob_dist)

#Part 1 Equations:
part1_equations = [
	"02x1+01x2->04x3 k=1",
	"01x1+02x3->03x2 k=2",
	"01x2+01x3->02x1 k=3"
]
#Exponentiation Reaction
exp_reactions = [
	"01q->01a1 k=1",
	"01a1+01y->01a1+02y' k=10000000000",
	"01a1->00o k=10000000",
	"01y'->01y k=1000"
]
#Logarithm Reactions
log2_reactions = [
	"01b->01a+01b k=1",
	"01a+02x->01c+01x'+01a k=10000000000",
	"02c->01c k=10000000000",
	"01a->00o k=10000000",
	"01x'->01x k=1000",
	"01c->01q k=1000"
]
#Problem 2 Reactions. ylog2(x)
p2_ylog2_x_reactions = log2_reactions + [
	"01q->01f k=100",
	"01y+01f->01f+01z+01y' k=10000000000",
	"01f->00o k=10000000",
	"01y'->01y k=1000"
]
p2_exp2_log2_x_reactions = exp_reactions + log2_reactions

p3_multiplication_reactions = abs_indicator("A1") + abs_indicator('B1') + abs_indicator('A2') + abs_indicator('B2') + [
	#Inputs: A1, A2, B1, B2
	#Split each input into 2:
	"01A1->01A1_0+01A1_1 K=" + FAST,
	"01B1->01B1_0+01B1_1 K=" + FAST,
	"01A2->01A2_0+01A2_1 K=" + FAST,
	"01B2->01B2_0+01B2_1 K=" + FAST,

	#Compute products C1, A1*B2, A2*B1, A2*B2
	#C1 = A1*B1: Takes A1_0 and B1_0 as input
	"01A1_ab+01B1_ab+01A1_0->01C1_mul k=" + SLOW,
	"01B1_0+01C1_mul->01C1_mul+01C1+01C1' k=" + VERY_FAST,
	"01C1_mul->00o k=" + FAST,
	"01C1'->01B1_0 k=" + MEDIUM,

	#A1*B2=D: Takes A1_1 and B2_0 as input
	"01A1_ab+01B2_ab+01A1_1->01D_mul k=" + SLOW,
	"01B2_0+01D_mul->01D_mul+01D+01D' k=" + VERY_FAST,
	"01D_mul->00o k=" + FAST,
	"01D'->01B2_0 k=" + MEDIUM,

	#A2*B2=E: Takes A2_1 and B2_1 as input
	"01A2_ab+01B2_ab+01A2_1->01E_mul k=" + SLOW,
	"01B2_1+01E_mul->01E_mul+01E+01E' k=" + VERY_FAST,
	"01E_mul->00o k=" + FAST,
	"01E'->01B2_1 k=" + MEDIUM,

	#A2*B1=F: Takes A2_0 and B1_1 as input
	"01A2_ab+01B1_ab+01A2_0->01F_mul k=" + SLOW,
	"01B1_1+01F_mul->01F_mul+01F+01F' k=" + VERY_FAST,
	"01F_mul->00o k=" + FAST,
	"01F'->01B1_1 k=" + MEDIUM,

	#Compute the addition reactions:
	"01D->01C2 k=" + VERY_FAST,
	"01E->01C2 k=" + VERY_FAST,
	"01F->01C2 k=" + VERY_FAST
]
p3_multiplication_counts = {
	'A1':5,
	'A2':10,
	'B1':15,
	'B2':20
}
p4_gcd_reactions = abs_indicator('x2') + abs_indicator('y2') + abs_indicator('x1') + abs_indicator('y1') + abs_indicator('x1') + abs_indicator('x3') + abs_indicator('y3') + abs_indicator('x') + abs_indicator('y') + [
	"01x3_ab+01y3_ab+01x->01x1+01x2 k=" + VERY_FAST,
	"01x3_ab+01y3_ab+01y->01y1+01y2 k=" + VERY_FAST,
	"01y_ab+01x_ab+01y2+01x2->00o k=" + MEDIUM,

	#If y>x:
	#"01x2_ab+01y1->00o k=" + FAST,
	#"01x2_ab+01y2->01y3 k=" + FAST
]
p4_gcd_counts = {
	'x':20,
	'y':50
}
p4_collatz_reactions = abs_indicator('e') + abs_indicator('O') + abs_indicator('x') + abs_indicator('x1') + abs_indicator('x3') + [
	"01e_ab+01O_ab+01x->01x1+01x2+01e_ab+01o_ab k=" + FAST, #Cloning Reaction
	#Even/Odd determining Reactions
	"01x_ab+02x2->01x_ab k=" + VERY_FAST, #Division by 2 until either one or no molecules remain. Used to determine evenness/oddness
	"01x_ab->01e k=" + MEDIUM, #Even indivator is produced in the absence of x
	"01e+01x2->01x2 k=" + VERY_FAST, #Even indicator is removed as long as x2 is still in the system
	"01x_ab+01x2->01x2+01O k=" + MEDIUM, #If there is a single X2 left, then create the odd indicator

	#Even odd detection is working.
	#Intermediate simulation results (even odd detection) (10000 iterations 10 times):
	#x = 30: ['x2: 0.0', "O_ab': 0.5", 'e: 687.0', "e_ab': 0.1", 'x1: 30.0', 'e_ab: 0.0', 'x_ab: 0.0', 'x: 0.0', 'O_ab: 1.1', 'O: 0.0', "x_ab': 0.5", 'o: 0.0', 'o_ab: 30.0']
	#x = 31: ['x2: 1.0', "x_ab': 0.8", "O_ab': 0.3", 'O_ab: 0.0', 'x: 0.0', 'x_ab: 0.0', 'O: 315.7', 'o: 0.0', "e_ab': 0.4", 'e: 0.0', 'o_ab: 31.0', 'e_ab: 1.1', 'x1: 31.0']

	#If the result is even, divide by 2
	"01e+02x1->01e+01x3 k=" + FAST, #Divide by two
	"01e+01x1_ab+01x3->01x+01e+01x1_ab k=" + FAST, #Once the division is done, react the x1 back into x

	#If the result is odd, multiple by 3 and add one:
	"01O+01x1->01O+03x3 k=" + FAST, #Multiply x1 by three and store it in x3
	"01O+01x1_ab+01x2->01x3+01O+01x1_ab k=" + VERY_FAST, #Change the last x2 into x3. The system should start producing e after this step.
	"01O+01x1_ab+01e+01x3->01x+01e+01O+01x1_ab k=" + FAST, #once x3 = 3x1 + 1 is finished, react the x3 back to x

	#Clean up and prepare for the next iteration
	"01x+01x3_abs+01e->00o k=" + VERY_FAST,
	"01x+01x3_abs+01O->00o k=" + VERY_FAST
]
p4_collatz_counts = {
	'x':3
}

if __name__ == "__main__":
	#Part 1:
	#p1_a_analyze_outcome(10000, 1000)
	p1_b_analyze_outcome(10000)

	#Part 2:
	#Z = y*log(x)
	#print(statistical_call_reaction(p2_ylog2_x_reactions, {'x':256, 'y':10, 'b':10}, 1000, 100))
	#print(statistical_call_reaction(p2_ylog2_x_reactions, {'x':16, 'y':100, 'b':10}, 1000, 100))
	#print(statistical_call_reaction(p2_ylog2_x_reactions, {'x':32, 'y':15, 'b':10}, 1000, 100))

	#Z = 2^log(x)
	#print(statistical_call_reaction(p2_exp2_log2_x_reactions, {'x':32,'y':1,'b':10}, 1000, 100))
	#print(statistical_call_reaction(p2_exp2_log2_x_reactions, {'x':64,'y':1,'b':10}, 1000, 100))
	#print(statistical_call_reaction(p2_exp2_log2_x_reactions, {'x':128,'y':1,'b':10}, 1000, 100))

	#Part 3:
	#print(statistical_call_reaction(p3_multiplication_reactions, p3_multiplication_counts, 10000, 10))

	#Part 4:
	#print(statistical_call_reaction(p4_collatz_reactions, p4_collatz_counts, 1000000, 1))
	#print(statistical_call_reaction(p4_gcd_reactions, p4_gcd_counts, 100000, 1))
