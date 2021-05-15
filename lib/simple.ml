open Printf
module L = Lazy
open Owl
open Defaults
open Biomechanics

let reuse = Cmdargs.check "-reuse"

(* read the setup files *)
let read_setup ~dynamic =
  let spontaneous, _, xstars, c =
    let file = if dynamic then "dynamic_setup.bin" else "setup.bin" in
    Miscellaneous.read_bin (Dir.in_dir file)
  in
  spontaneous, xstars, c


let in_dir s =
  let prefix = Cmdargs.(get_string "-prefix" |> default "") in
  Dir.in_dir Printf.(sprintf "%s/%s" prefix s)


(* both spontaneous and c are identical for setup and dynamic setup *)
let spontaneous, _, c = read_setup ~dynamic:false

module Prms = struct
  include Defaults.Main

  let spontaneous = spontaneous
end

include Prms
module D = Dynamics.Make (Prms)

let q_mov c =
  let module G =
    Gramians.Make (struct
      let a = Defaults.Main.a
      let m1_slice = Prms.m1_slice
      let gamma = Some Mat.(c *@ Defaults.Main.gamma)
    end)
  in
  let q = Miscellaneous.unit_trace G.OSub.m_norm in
  (* symmetrize for numerical stability *)
  Mat.(0.5 $* q + transpose q)


let q_prep c =
  let cg = Mat.(c *@ gamma) in
  Miscellaneous.unit_trace Mat.(transpose cg *@ cg)


let save ?(gate_torque_off = false) ~qs ~xstar (module P : Tspecs.T) ~path (x, y, z) =
  let open P in
  Mat.(save_txt (transpose x) ~out:(path "x"));
  Mat.(save_txt (transpose (x - spontaneous)) ~out:(path "x_relative"));
  Miscellaneous.maybe y (fun y -> Mat.(save_txt (transpose y) ~out:(path "y"))) ();
  Miscellaneous.maybe z (fun z -> Mat.(save_txt (transpose z) ~out:(path "z"))) ();
  (* gate the torques off until movement onset *)
  let torques = Mat.(c *@ gamma *@ x) in
  let torque_gate_on = Owl.Maths.(round ((off1 +. on) /. sampling_dt)) |> int_of_float in
  (if (not gate_torque_off) && torque_gate_on > 0
  then
    Mat.(
      torques.${[]; [ 0; pred torque_gate_on ]} <- zeros (row_num torques) torque_gate_on));
  let torques = Mat.transpose torques in
  Mat.(save_txt torques ~out:(path "torques"));
  (* compute arm trajectory *)
  let theta = Arm.theta_trajectory ~dt:sampling_dt torques in
  Mat.(save_txt (Arm.unpack_sequence theta) ~out:(path "theta"));
  let hand = Array.map Arm.hand_of theta in
  Mat.(save_txt (Arm.unpack_sequence hand) ~out:(path "hand"));
  (* compute the timecourse of the quadratic error *)
  let cost q = Mat.(transpose (x - xstar) *@ q *@ (x - xstar) |> diag |> transpose) in
  Array.map cost qs |> Mat.concatenate ~axis:1 |> Mat.save_txt ~out:(path "cost")


(* --------------------------------------------------------------------------------
   --     VANILLA LQR FEEDBACK                                                 --
   -------------------------------------------------------------------------------- *)

let vanilla
    ?(verbose = true)
    ?(reuse = reuse)
    ?(r2 = r2_vanilla)
    ?(x0 = Some spontaneous)
    ?nl
    ~q
    ~xstars
    (module P : Tspecs.T)
  =
  let open P in
  (* compute the optimal feedback matrix *)
  let ((xy, yx) as lqr_feedback) =
    Miscellaneous.maybe_reuse
      reuse
      "vanilla_feedback.bin"
      (lazy (Controller.vanilla ~a ~pmd_slice ~r2 ~q))
  in
  (* simulate the dynamics for each xstar *)
  ( Array.mapi
      (fun j xstars ->
        Array.mapi
          (fun i xstar ->
            if verbose then printf "Vanilla LOOP r%i %i\n%!" (succ j) (i + 1);
            let x =
              D.simulate_instantaneous_feedback
                ~sampling_dt
                ~duration
                ~mov_input
                ~x0
                ~verbose
                ?nl
                ~loop:(lqr_feedback, gain_f)
                xstar
            in
            x, None, None)
          xstars)
      xstars
  , Some (xy, yx, None) )


(* --------------------------------------------------------------------------------
   --     DYNAMIC FEEDBACK                                                       --
   -------------------------------------------------------------------------------- *)

let dynamic ?(verbose = true) ?(x0 = Some spontaneous) ~xstars (module P : Tspecs.T) =
  let open P in
  (* compute the optimal feedback matrix *)
  let _, ((xz, zy, yx) as decomposed) =
    Miscellaneous.read_bin (Dir.in_dir "dynamic_feedback.bin")
  in
  (* simulate the dynamics for each xstar *)
  ( Array.mapi
      (fun j xstars ->
        Array.mapi
          (fun i xstar ->
            if verbose then printf "DYNAMIC LOOP r%i %i\n%!" (succ j) (i + 1);
            let x, y, z =
              D.simulate_dynamic_feedback
                ~sampling_dt
                ~duration
                ~mov_input
                ~x0
                ~verbose
                ~taus
                ~loop:(decomposed, gain_f, input_gain_f)
                xstar
            in
            x, Some y, Some z)
          xstars)
      xstars
  , Some (xz, zy, Some yx) )


(* --------------------------------------------------------------------------------
   --     NAIVE CONTROL STRATEGY WITH SPATIALLY CONSTANT INPUTS                  --
   -------------------------------------------------------------------------------- *)

let naive ?(verbose = true) ?(x0 = Some spontaneous) ~xstars (module P : Tspecs.T) =
  let open P in
  (* zero feedback gain *)
  let (_ as no_feedback) = Mat.zeros 1 n, Mat.zeros n 1 in
  (* simulate the dynamics for each xstar *)
  ( Array.mapi
      (fun j xstars ->
        Array.mapi
          (fun i xstar ->
            if verbose then printf "Native LOOP r%i %i\n%!" (succ j) (i + 1);
            let x =
              D.simulate_instantaneous_feedback
                ~sampling_dt
                ~duration
                ~mov_input
                ~x0
                ~verbose
                ~loop:(no_feedback, gain_f)
                xstar
            in
            x, None, None)
          xstars)
      xstars
  , None )
