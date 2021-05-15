open Owl

module type PT = sig
  val n : int
  val n_e : int
  val n_i : int
  val w_rec : Mat.mat
  val spontaneous : Mat.mat
  val dt : float
  val tau : float (* time constant of cortical dynamics *)

  val nl : Mat.mat -> Mat.mat
  val noise_prms : (float * Mat.mat * (float -> float)) option

  (* optional input noise parameters *)
end

module Make (P : PT) : sig
  include PT

  val a : Mat.mat

  (** steady input *)
  val h : Mat.mat

  (* you can run multiple trials in parallel, by specifying an [n x n_trials]
     array of initial conditions for x0;
     the result is either a neuron x time matrix if n_trials = 1,
     or an n_trials x neuron x time ndarray otherwise *)
  val simulate_instantaneous_feedback
    :  ?eta0:Arr.arr
    -> sampling_dt:float
    -> duration:float
    -> mov_input:(float -> float)
    -> x0:Mat.mat option
    -> ?verbose:bool
    -> ?perturbation:(float -> Mat.mat)
    -> ?nl:(Mat.mat -> Mat.mat)
    -> loop:(Mat.mat * Mat.mat) * (float -> float)
    -> (* loop and temporal loop gain *)
       Mat.mat
    -> (* xstar *)
       Arr.arr

  val simulate_dynamic_feedback
    :  ?eta0:Arr.arr
    -> sampling_dt:float
    -> duration:float
    -> mov_input:(float -> float)
    -> x0:Mat.mat option
    -> ?verbose:bool
    -> ?perturbation:(float -> Mat.mat)
    -> ?nl:(Mat.mat -> Mat.mat)
    -> taus:float * float
    -> loop:(Mat.mat * Mat.mat * Mat.mat) * (float -> float) * (float -> float)
    -> (* loop and temporal loop gain *)
       Mat.mat
    -> (* xstar *)
       Arr.arr * Arr.arr * Arr.arr
end
