open Owl

(* read the setup files *)
val read_setup : dynamic:bool -> Mat.mat * Mat.mat array array * Mat.mat

(* save with prefix *)
val in_dir : string -> string

(* some default stuff *)
val c : Mat.mat
val spontaneous : Mat.mat

(* relevant metrics *)
val q_mov : Mat.mat -> Mat.mat
val q_prep : Mat.mat -> Mat.mat

val save
  :  ?gate_torque_off:bool
  -> qs:Mat.mat array
  -> xstar:Mat.mat
  -> (module Tspecs.T)
  -> path:(string -> string)
  -> Arr.arr * Arr.arr option * Arr.arr option
  -> unit

val vanilla
  :  ?verbose:bool
  -> ?reuse:bool
  -> ?r2:float
  -> ?x0:Mat.mat option
  -> ?nl:(Mat.mat -> Mat.mat)
  -> q:Mat.mat
  -> xstars:Mat.mat array array
  -> (module Tspecs.T)
  -> (Arr.arr * Arr.arr option * Arr.arr option) array array
     * (Mat.mat * Mat.mat * Mat.mat option) option

val dynamic
  :  ?verbose:bool
  -> ?x0:Owl.Mat.mat option
  -> xstars:Mat.mat array array
  -> (module Tspecs.T)
  -> (Arr.arr * Arr.arr option * Arr.arr option) array array
     * (Mat.mat * Mat.mat * Mat.mat option) option

val naive
  :  ?verbose:bool
  -> ?x0:Mat.mat option
  -> xstars:Mat.mat array array
  -> (module Tspecs.T)
  -> (Arr.arr * Arr.arr option * Arr.arr option) array array
     * (Mat.mat * Mat.mat * Mat.mat option) option
