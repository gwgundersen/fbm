/* MC.H - Interface to Markov chain Monte Carlo module. */

/* Copyright (c) 1995, 1996, 1998 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, or modify this program 
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
 */


/* OPERATIONS TO PERFORM EACH ITERATION.  This array of structures lists
   the operations to be performed (in sequence) in each iteration of a
   Markov chain Monte Carlo procedure.  The type of the operation is specified
   by a one-character identifier, with zero meaning that the slot is empty.  
   The other fields are meaningful only when the type is appropriate. 
   Some operations apply to a group of other operatons, in which case the
   end of the group is marked by an E (or by the end of the whole list). 

   Stored in log files under type 'o'.  Changes may invalidate old log files. */

#define Max_mc_ops 20	/* Maximum number of operations in an iteration */

#define Group_ops "Rt"	/* Operations that are followed by a group of ops */

typedef struct 
{ 
  struct 		/* List of operations to perform */
  { 
    int type;		  /* Type of operation, zero for empty slot */
    int repeat_count;	  /* Repetition count for type 'R'' */

    int steps;		  /* Number of steps in trajectory, or max steps, or
                             max intervals for slice sampling */

    float stepsize_adjust;/* Adjustment factor for stepsizes */
    float stepsize_alpha; /* Gamma param for stepsize dist, zero is infinity */

    int window; 	  /* Window size for hybrid Monte Carlo updates */
    int jump;	  	  /* Steps in each jump for hybrid Monte Carlo */

    float heatbath_decay; /* Momentum decay for heatbath step */

    float temper_factor;  /* Tempering factor for tempered hybrid Monte Carlo */
    float app_param;	  /* Parameter for application-specific procedure */

    int firsti, lasti;	  /* Indexes of first and last coordinates to apply
                             single-variable updates to */

    float refresh_prob;	  /* Prob. of refresh in overrelaxed slice sampling */

    int refinements;	  /* Number of refinement steps */

    int in_steps;	  /* Maximum number of inside steps in trajectory,
			     zero if this feature is not being used */

    float app_param2;	  /* Second application-specific parameter */

    int reserved[2];      /* Reserved for future use */

    char appl[101];	  /* Name of application-specific procedure */

  } op[Max_mc_ops];

} mc_ops;


/* METHOD FOR COMPUTING TRAJECTORIES.  This structure describes how dynamical
   trajectories should be computed. 

   Stored in log files under type 't'.  Changes may invalidate old log files. */

#define Max_approx 100	/* Maximum number of energy approximations */

typedef struct
{ 
  int type;		/* Type of discretization, currently always 'L' */
  int halfp;		/* Should first and last half steps be for p? */
  int N_approx;		/* Number of approximations to the energy function;
			     when negative, each is used twice, symmetrically */

  int reserved[5];	/* Reserved for future use */

} mc_traj;


/* TEMPERING SCHEDULE AND BIASES.  This structure gives the list of 
   inverse temperatures used in tempering schemes, together with the
   biases associated with each temperature used in simulated tempering.

   Stored in log files under type 'm'.  Changes may invalidate old log files. 

   Note that the "temperatures" here are different from the "temperatures"
   used in connection with heatbath updates, and for acceptance of
   proposed moves.  One might fiddle with both temperatures simultaneously,
   though this would perhaps be a bit silly; normally at least one of these
   temperatures would be fixed at one. */

#define Max_temps 1001	/* Maximum number of temperatures in a schedule */

typedef struct
{
  struct mc_temp_bias	/* Tempering schedule, ending with inv_temp==1 */
  { 
    float inv_temp;	  /* Inverse temperature */
    float bias;		  /* Log of bias factor */

  } sched[Max_temps];

} mc_temp_sched;


/* TEMPERING STATE.  This structure contains the current inverse temperature
   and the direction in which the temperature is to be changed (not always
   relevant). 

   Stored in log files under type 'b'.  Changes may invalidate old log files. */

typedef struct
{
  float inv_temp;	/* Current inverse temperature, must correspond to an 
			     entry in the tempering schedule, except that it's
                             set to zero at the start of an ais run. */

  int temp_dir;		/* Direction in which to change inverse temperature   */

} mc_temp_state;


/* INFO ON MONTE CARLO ITERATION.  This structure records various bits of 
   information concerning the current iteration.  The temperature and decay
   values are derived from user specifications; the approx_order field is
   part of the simulation state; the stepsize_factor field is selected at
   random each iteration; the remaining fields reflect the results. 

   Stored in log files under type 'i'.  Changes may invalidate old log files. */

typedef struct
{
  float temperature;	/* Temperature used during this iteration */
  float decay;		/* Heatbath decay used during this iteration */

  char approx_order[Max_approx]; /* Order in which approximations are used */

  float stepsize_factor;/* Factor last used to adjust stepsizes */

  double delta;		/* Change in total energy for last state proposed */
  int move_point;	/* Last point moved to along trajectory */
  int window_offset;	/* Offset of start state within window, for hybrid */
  int rejects;		/* Number of rejections in this iteration */
  int proposals;	/* Number of proposals in this iteration */
  int slice_calls;	/* Number of calls of slice sampling procedures */
  int slice_evals;	/* Number of energy evaluations in slice calls */
  int time;             /* Cumulative cpu-usage in ms */

  float log_weight;	/* Log of weight for importance sampling */
  float log_tt_weight;	/* Log of weight from last tempered transition */
  float log_tt_weight2;	/* Log of combined weight from last tempered trans. */

  int reserved[1];	/* Reserved for future use */

} mc_iter;


/* DYNAMICAL VARIABLES.  This structure contains pointers to the position
   and momentum components of the state, as well to the stepsizes to use
   for each component and saved values of the potential and kinetic energy 
   (useful to avoid unnecessary recomputation).  The inverse temperature
   used in simulated tempering is also stored here. */

typedef double mc_value;/* Type of position and momentum values */

typedef struct
{ 
  int dim;		/* Dimensionality of position and momentum */
  int aux_dim;		/* Dimensionality of auxiliary variables */

  mc_value *q;		/* Position variables */
  mc_value *p;		/* Momentum variables */
  mc_value *aux;	/* Auxiliary variables */

  mc_temp_state *temp_state; /* State for simulated tempering */
  int temp_index;	/* Index of inverse temperature in schedule */

  mc_value *grad;	/* Gradient of potential energy w.r.t. position */
  mc_value *stepsize;	/* Stepsizes to use for each component */

  double pot_energy;	/* Potential energy of position */
  double kinetic_energy;/* Kinetic energy of momentum */

  int know_pot;		/* Is potential energy up-to-date? */
  int know_kinetic;	/* Is kinetic energy up-to-date? */

  int know_grad;	/* Which gradient approx. is known?  Zero if none. */

} mc_dynamic_state;


/* PROCEDURES PROVIDED BY THE APPLICATION. */

extern void mc_app_record_sizes (log_gobbled *);
extern void mc_app_initialize (log_gobbled *, mc_dynamic_state *);
extern void mc_app_save (mc_dynamic_state *, log_file *, int);

extern int mc_app_sample (mc_dynamic_state *, char *, double, double,
                          mc_iter *, mc_temp_sched *);

extern void mc_app_energy (mc_dynamic_state *, int, int, double *, mc_value *);
extern int mc_app_zero_gen (mc_dynamic_state *);

extern void mc_app_stepsizes (mc_dynamic_state *);


/* MARKOV CHAIN MONTE CARLO PROCEDURES. */

void mc_iter_init  (mc_dynamic_state *, mc_ops *, mc_traj *, mc_temp_sched *);
void mc_iteration  (mc_dynamic_state *, mc_iter *, log_gobbled *, void *, int);

void mc_traj_init    (mc_traj *, mc_iter *);
void mc_traj_permute (void);
void mc_trajectory   (mc_dynamic_state *, int, int);

void mc_record_sizes     (log_gobbled *);
void mc_heatbath         (mc_dynamic_state *, float, float);
void mc_radial_heatbath  (mc_dynamic_state *, float);
double mc_kinetic_energy (mc_dynamic_state *);
void mc_temp_present     (mc_dynamic_state *, mc_temp_sched *);
int mc_temp_index        (mc_temp_sched *, float);
double mc_energy_diff    (mc_dynamic_state *, mc_temp_sched *, int);

void mc_value_copy (mc_value *, mc_value *, int);
