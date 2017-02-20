#Molecular Logic Elements for EE 5393 HW 1
#Uses reaction_simulation.py
#Author: Owen Hoffend

from reaction_simulation import parse_reactions
from reaction_simulation import stochastic_sim

def mux(iterations):
	reactions = [
		"01x+01S1->01S1+01z k=1",
		"01y+01S2->01S2+01z k=1"
	]
	counts = {
		'x': 50,
		'y': 10,
		'S1': 10
	}
	sim_input = parse_reactions(reactions, counts)
	print(sim_input)
	output = stochastic_sim(sim_input[1], sim_input[0], iterations)
	return output

def no_products(iterations):
	reactions = ['01x->00y k=1']
	counts = {'x': 100}

	sim_input = parse_reactions(reactions, counts)
	print(sim_input)
	output = stochastic_sim(sim_input[1], sim_input[0], iterations)
	return output

if __name__ == "__main__":
	print(no_products(100))
