(* dynamics of the movement phase, written under Algodiff
   to be used in the optimization of C and xstars *)

open Owl

module type PT = sig
  val dt : float
  val sampling_dt : float
  val tau : float
  val spontaneous : Mat.mat
  val mov_input : float -> float
end

module Make (P : PT) : sig
  val run
    :  ?nl:(Algodiff.D.t -> Algodiff.D.t)
    -> w_rec:Algodiff.D.t
    -> n_bins:int
    -> layers:((float * Mat.mat * Algodiff.D.t) * (float * Mat.mat * Algodiff.D.t)) option
    -> Algodiff.D.t
    -> Algodiff.D.t array
end
