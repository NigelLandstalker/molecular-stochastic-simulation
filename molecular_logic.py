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
	return averages

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
	"01x->01a k=1",
	"01a+01y->01a+02y' k=10000000000",
	"01a->00o k=10000000",
	"01y'->01y k=1000"
]
exp_counts = {
	'x':5,
	'y':1
}

#Logarithm Reactions
log2_reactions = [
	"01b->01a+01b k=1",
	"01a+02x->01c+01x'+01a k=10000000000",
	"02c->01c k=10000000000",
	"01a->00o k=10000000",
	"01x'->01x k=1000",
	"01c->01y k=1000"
]
log2_counts = {
	'b':10,
	'x':128
}

if __name__ == "__main__":
	print(statistical_call_reaction(log2_reactions,log2_counts, 10000, 10))
