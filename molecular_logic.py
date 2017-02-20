#Molecular Logic Elements for EE 5393 HW 1
#Uses reaction_simulation.py
#Author: Owen Hoffend

from reaction_simulation import parse_reactions
from reaction_simulation import stochastic_sim

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
	new_output = [key + str(outputs[unique_molecules[key]]) for key in unique_molecules]
	return new_output

#Mux Reaction (as a test)
mux_reactions = [
	"01x+01S1->01S1+01z k=1",
	"01y+01S2->01S2+01z k=1"
]
mux_counts = {
	'x': 50,
	'y': 10,
	'S1': 10
}

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
	#Absence Indicator Reactions:
	#"00o->01e k=1",
	#"01x+01e->01x k=100",
	#"02e->01e k=100",
	#Multiplication Reactions
	#"01q+01e->01f k=100",
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

if __name__ == "__main__":
	print(statistical_call_reaction(p2_exp2_log2_x_reactions, p2_exp2_log2_x_counts, 100000, 10))
	#print(statistical_call_reaction(p2_ylog2_x_reactions, p2_ylog2_x_counts, 1000, 100))
