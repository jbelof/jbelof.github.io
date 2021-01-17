/************************************************************************/
/*									*/
/* This code implements the stochastic Markov model of Astumian		*/
/* and it's (seemingly) paradoxical behaviour.  Games 0 are 1 are both	*/
/* both probalistically losing, yet randomly switching between the two	*/
/* games (denoted Game 2) yields a winning outcome on average.		*/
/* The initial condition is the middle state, i=2 (phi=0.5)		*/
/*									*/
/* An interesting addition to this story is that using a		*/
/* correlated/biased fluctuation spectrum (denoted Game 3) results in	*/
/* net losing in excess of simply playing a single game!		*/
/*									*/
/* We describe each state point by an order parameter, which is either	*/
/* 0, 0.25, 0.5, 0.75 or 1.  When drawing a random number to decide	*/
/* which game to play, rather than draw a RN from a uniform		*/
/* distribution, [0,1], we draw from [0,phi]				*/
/*									*/
/* D.R. Astumian, Am. J. Phys. 2005, 73(2):178-183			*/
/*									*/
/* compilation:								*/
/*	gcc -o a astumian.c -lm						*/
/*									*/
/* usage:								*/
/*	./a: <game> <numsteps>						*/
/*		<game> - integer game-type of 0, 1,			*/
/*			2 (uniform) or 3 (correlated)			*/
/*		<numsteps> - integer number of steps to perform		*/
/*									*/
/* @2009, Jonathan Belof						*/
/*									*/
/************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define SEED		2358

#define LOSE_STATE	0
#define WIN_STATE	4
#define INITIAL_STATE	2

#define NUM_GAMES	2
#define NUM_STATES	5


/* initialize the PRNG */
void seed_rng(long int seed) {

	srand48(seed);

}

/* uses standard UNIX 48-bit linear congruence PSRG */
double get_random_number(void) {

	return(drand48());

}

void usage(char *progname) {

	fprintf(stderr, "%s: <game> <numsteps> \n", progname);
	fprintf(stderr, "\t<game> - integer game-type of 0,1, 2 (uniform) or 3 (correlated)\n");
	fprintf(stderr, "\t<numsteps> - integer number of steps to perform\n");
	exit(1);

}


int main(int argc, char **argv) {

	int game, numsteps;					/* simulation parameters */
	int uniform = 0, correlated = 0;			/* random switching flags */
	int state, losses, wins;				/* game state and outcome */
	int step;						/* loop index */
	double game_mm[NUM_GAMES][NUM_STATES][NUM_STATES];	/* Markov transition matrix */
	int nextstate, trans;					/* the next state to transition to, transition under consideration */
	double cummulant;					/* working variable for trans. prob. interval calc. */
	double rand_trans;					/* random number for transition probability */
	double phi;						/* order parameter, 0=loss and 1=win */

	/* read arguments */
	if(argc != 3) usage(argv[0]);

	game = atoi(argv[1]);
	numsteps = atoi(argv[2]);

	/* argument checking */
	if ((game < 0) || (game > 3)) {
		fprintf(stderr, "error: invalid game-type requested, should be 0, 1, 2 or 3\n");
		usage(argv[0]);
	} else printf("Running game-type %d\n", game);

	/* these flags specify random switching between games */
	if(game == 2) {
		game = 0;
		uniform = 1;
	} else if(game == 3) {
		game = 0;
		correlated = 1;
	}

	if(numsteps < 0) {
		fprintf(stderr, "error: invalid number of simulation steps provided\n");
		usage(argv[0]);
	} else printf("Running %d simulation steps\n", numsteps);


	/* seed RNG */
	seed_rng(SEED);


	/* setup the Markov matrix for game 0 */
	game_mm[0][0][0] = 0;
	game_mm[0][0][1] = 0;
	game_mm[0][0][2] = 0;
	game_mm[0][0][3] = 0;
	game_mm[0][0][4] = 0;
	game_mm[0][1][0] = 4./36.;
	game_mm[0][1][1] = 24./36.;
	game_mm[0][1][2] = 8./36.;
	game_mm[0][1][3] = 0;
	game_mm[0][1][4] = 0;
	game_mm[0][2][0] = 0;
	game_mm[0][2][1] = 5./36.;
	game_mm[0][2][2] = 29./36.;
	game_mm[0][2][3] = 2./36.;
	game_mm[0][2][4] = 0;
	game_mm[0][3][0] = 0;
	game_mm[0][3][1] = 0;
	game_mm[0][3][2] = 4./36.;
	game_mm[0][3][3] = 24./36.;
	game_mm[0][3][4] = 8./36.;
	game_mm[0][4][0] = 0;
	game_mm[0][4][1] = 0;
	game_mm[0][4][2] = 0;
	game_mm[0][4][3] = 0;
	game_mm[0][4][4] = 0;


	/* setup the Markov matrix for game 1 */
	game_mm[1][0][0] = 0;
	game_mm[1][0][1] = 0;
	game_mm[1][0][2] = 0;
	game_mm[1][0][3] = 0;
	game_mm[1][0][4] = 0;
	game_mm[1][1][0] = 5./36.;
	game_mm[1][1][1] = 29./36.;
	game_mm[1][1][2] = 2./36.;
	game_mm[1][1][3] = 0;
	game_mm[1][1][4] = 0;
	game_mm[1][2][0] = 0;
	game_mm[1][2][1] = 4./36.;
	game_mm[1][2][2] = 24./36.;
	game_mm[1][2][3] = 8./36.;
	game_mm[1][2][4] = 0;
	game_mm[1][3][0] = 0;
	game_mm[1][3][1] = 0;
	game_mm[1][3][2] = 5./36.;
	game_mm[1][3][3] = 29./36.;
	game_mm[1][3][4] = 2./36.;
	game_mm[1][4][0] = 0;
	game_mm[1][4][1] = 0;
	game_mm[1][4][2] = 0;
	game_mm[1][4][3] = 0;
	game_mm[1][4][4] = 0;

	/* initialize the system */
	losses = 0, wins = 0;

	/* the main loop of the simulation */
	for(step = 0; step < numsteps; step++) { 	/* loop over each game */

		/* keep running the stochastic walk until we get a result */
		state = INITIAL_STATE;
		while(1) {

			/* determine the next state in the Markov chain */
			rand_trans = get_random_number();
			for(trans = 0, nextstate = -1, cummulant = 0; trans < NUM_STATES; trans++) {

				cummulant += game_mm[game][state][trans];
				if(nextstate < 0) {
					if(rand_trans < cummulant) {
						nextstate = trans;
						break;
					}
				}

			}

			/* now let's transition to the new state */
			state = nextstate;

			/* randomly switch MTM */
			if(uniform) {

				game = (int)rint(get_random_number());

			} else if (correlated) {

				phi = ((double)state)/((double)(NUM_STATES-1));
				game = (int)rint(phi*get_random_number());

			}

			/* if we reached the end, score our result and reset the walk */
			if(state == LOSE_STATE) {
				++losses; break;
			} else if (state == WIN_STATE) {
				++wins; break;
			}

		}

	}

	/* output the results */
	printf("losses / wins = %d/%d = %f\n", losses, wins, ((double)losses)/((double)wins));
	printf("analytic: game0 = %f, game1 = %f, uniform = %f\n", (20./16.), (20./16.), (81./100.));


	exit(0);

}


