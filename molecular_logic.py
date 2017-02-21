#Molecular Logic Elements for EE 5393 HW 1
#Uses reaction_simulation.py
#Author: Owen Hoffend

from reaction_simulation import parse_reactions
from reaction_simulation import stochastic_sim

SLOW = '1'
MEDIUM = '100'
FAST = '10000'
VERY_FAST = '1000000'

def call_reaction(reaction, counts, iterations):
	sim_input = parse_reactions(reaction, counts)
	output = stochastic_sim(sim_input[1], sim_input[0], iterations)
	return output

def statistical_call_reaction(reaction, counts, iterations, trials):
	sim_input = parse_reactions(reaction, counts)
	outputs = [stochastic_sim(sim_input[1], sim_input[0], iterations) for _ in range(trials)]
	averages = [(sum(counts) / trials) for counts in list(zip(*outputs))]
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

#Exponentiation Reaction
exp_reactions = [
	"01q->01a1 k=1",
	"01a1+01y->01a1+02y' k=10000000000",
	"01a1->00o k=10000000",
	"01y'->01y k=1000"
]
exp_counts = {
	'q':5,
	'y':1
}

#Logarithm Reactions
log2_reactions = [
	"01b->01a+01b k=1",
	"01a+02x->01c+01x'+01a k=10000000000",
	"02c->01c k=10000000000",
	"01a->00o k=10000000",
	"01x'->01x k=1000",
	"01c->01q k=1000"
]
log2_counts = {
	'b':10,
	'x':128
}
#Problem 2 Reactions. ylog2(x)
p2_ylog2_x_reactions = log2_reactions + [
	"01q->01f k=100",
	"01y+01f->01f+01z+01y' k=10000000000",
	"01f->00o k=10000000",
	"01y'->01y k=1000"
]
p2_ylog2_x_counts = {
	'b':10,
	'x':256, #Value to modify
	'y':10 #Value to modify
}

p2_exp2_log2_x_reactions = exp_reactions[1:] + log2_reactions + [
	#Modified exp input: waits for absence indicator to fire.
	"01q+01e->01a1 k=1",
	#Absence Indicator Reactions
	"00o->01e k=1",
	"01x+01e->01x k=100",
	"02e->01e k=100"
]
p2_exp2_log2_x_counts = {
	'y':1,
	'b':10,
	'x':32 #Value to modify
}

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

p4_gcd_reactions = [

]
p4_gcd_counts = {

}

p4_collatz_reactions = [

]
p4_collatz_counts = {

}

if __name__ == "__main__":
	#print(statistical_call_reaction(addition_reactions,addition_counts, 1000, 100))
	#print(statistical_call_reaction(p2_ylog2_x_reactions, p2_ylog2_x_counts, 1000, 100))
	print(statistical_call_reaction(p3_multiplication_reactions, p3_multiplication_counts, 10000, 10))
	#print(abs_indicator('q'))
