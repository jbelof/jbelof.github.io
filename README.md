# Astumian's Game

This code implements the stochastic Markov model of Astumian [D.R. Astumian, Am. J. Phys. 2005, 73(2):178-183] and it's (seemingly) paradoxical behaviour.  Games labelled as 0 or 1 are both both probalistically losing, yet randomly switching between the two games (denoted Game 2) yields a winning outcome on average. The initial condition is the middle state, i=2 (phi=0.5) and transitions
proceed until either losing (0/phi=0) or winning (4/phi=1) states are reached.

An interesting addition to this story is that using a correlated/biased fluctuation spectrum (denoted Game 3) results in net losing in excess of simply playing a single game!

We describe each state point by an order parameter, which is either 0, 0.25, 0.5, 0.75 or 1.  When drawing a random number to decide which game to play, rather than draw a RN from a uniform distribution, [0,1], we draw from [0,phi]

The results of this toy model highlight how non-equilibrium systems, subject to thermal noise (uncorrelated or correlated) driving a system across a rough energy landscape, can drive reactions in reverse. 

@2009, Jonathan Belof

## Getting Started

After obtaining the source code, please consult the Makefile to set any specific compiler flags (defaults are gcc with std C library).

The code implements direct numerical simulation of a Markov Chain Monte Carlo for two different reaction networks, with specific barriers and energy wells, such that "playing" either one of them leads to a "losing" game on average. Randomly alternating between games leads to a "winning" result.


## Installing

Compilation is simple and relies on only standard libraries:

$ make
gcc -c -O3 -DDEBUG -I. astumian_game.c
gcc -O3 -DDEBUG *.o -o ag


## Running the examples

Run the binary without arguments to obtain the usage of command line input:

$ ./ag
./ag: <game> <numsteps> 
	<game> - integer game-type of 0,1, 2 (uniform) or 3 (correlated)
	<numsteps> - integer number of steps to perform

To run the (losing) Game 0:

$ ./ag 0 1000000
Running game-type 0
Running 1000000 simulation steps
losses / wins = 555560/444440 = 1.250023
analytic: game0 = 1.250000, game1 = 1.250000, uniform = 0.810000

where the MCMC results are shown in the penultimate line, compared with the clsed form analytic result obtain via product of state probabilities.

Likewise for running Game 1:

$ ./ag 1 1000000
Running game-type 1
Running 1000000 simulation steps
losses / wins = 555588/444412 = 1.250164
analytic: game0 = 1.250000, game1 = 1.250000, uniform = 0.810000

with the explicit numerical calculation approaching the analytic result with more instantiations.

With Game 2, we flip a coin and randomly select Game 1 or Game 2 Markov probabilities at each transition point; playing this hybride game leads to a reversal of the system!

$ ./ag 2 1000000
Running game-type 2
Running 1000000 simulation steps
losses / wins = 446432/553568 = 0.806463
analytic: game0 = 1.250000, game1 = 1.250000, uniform = 0.810000

(note that the analytic result to compare with is the 3rd one called "uniform").


## Authors

* **Jon Belof** [jbelof@github](https://github.com/jbelof)  

[home page](http://people.llnl.gov/belof1)
[google scholar](https://scholar.google.com/citations?user=gNrlNbwAAAAJ&hl=en)  
[research gate](https://www.researchgate.net/profile/Jon_Belof)  
[linkedin](http://www.linkedin.com/in/jbelof)  
[web profile](http://jbelof.academia.edu)  


## License

This project is licensed under the GNU General Public License v3, please see GPL_license.txt for details.


