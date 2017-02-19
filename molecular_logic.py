#Molecular Logic Elements for EE 5393 HW 1
#Uses reaction_simulation.py
#Author: Owen Hoffend

import reaction_simulation

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
	sim_input = reaction_simulation.parse_reactions(reactions, counts)
	print(sim_input)
	output = reaction_simulation.stochastic_sim(sim_input[1], sim_input[0], iterations)
	return output

if __name__ == "__main__":
	print(mux(1000))
