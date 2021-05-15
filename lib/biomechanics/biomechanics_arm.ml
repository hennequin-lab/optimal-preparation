open Owl

module Prms = struct
  type 'a t = { torques : 'a } [@@deriving prms]
end

open Prms
module O = Owl_opt_lbfgs.D.Make (Prms)

type position
type angular

(* angular or hand state *)
type 'a state =
  { x1 : Algodiff.D.t
  ; x1_dot : Algodiff.D.t
  ; x2 : Algodiff.D.t
  ; x2_dot : Algodiff.D.t
  }

(* ----------------------------------------------------------------------------
   ---     SOME UTILITY FUNCTIONS                                           ---
   ---------------------------------------------------------------------------- *)

let unpack x = Algodiff.D.(unpack_flt (primal x))

let unpack_state state =
  Mat.of_arrays
    [| [| unpack state.x1; unpack state.x1_dot; unpack state.x2; unpack state.x2_dot |] |]


let unpack_sequence states = Array.map unpack_state states |> Mat.concatenate ~axis:0
let deg x = x *. Const.pi /. 180.

(* ----------------------------------------------------------------------------
   ---     CONSTANT PARAMETERS                                              ---
   ---------------------------------------------------------------------------- *)

(* arm lengths *)
let _L1 = 0.3
let _L2 = 0.33

(* moments of inertia *)
let _I1 = 0.025
let _I2 = 0.045

(* masses *)
let _M1 = 1.4
let _M2 = 1.0

(* B matrix *)
let _B11 = 0.05
let _B12 = 0.025
let _B21 = 0.025
let _B22 = 0.05

(* distances *)
let _S1 = 0.11
let _S2 = 0.16

(* initial theta *)
let _THETA0 =
  (* let theta1 = 20. *. Const.pi /. 180. in  (* for 27 conditions *) *)
  let theta1 = 10. *. Const.pi /. 180. in
  (* for 8 conditions *)
  let joint_x = _L1 *. cos theta1 in
  (* we want L2 * cos (theta1 + theta2) = - joint_x *)
  let theta2 = acos (-.joint_x /. _L2) -. theta1 in
  Algodiff.D.{ x1 = F theta1; x1_dot = F 0.; x2 = F theta2; x2_dot = F 0. }


(* derived constants *)
let _A1 = _I1 +. _I2 +. (_M2 *. Maths.sqr _L1)
let _A2 = _M2 *. _L1 *. _S2
let _A3 = _I2

(* ----------------------------------------------------------------------------
   ---     KINEMATICS MATH FUNCTIONS                                        ---
   ---------------------------------------------------------------------------- *)

(* forward dynamics *)
let theta_dot_dot =
  let open Algodiff.D in
  let b = Maths.of_arrays [| [| F _B11; F _B12 |]; [| F _B21; F _B22 |] |] in
  fun theta torque ->
    let torque = Maths.transpose torque in
    let m =
      let z = Maths.(F _A2 * cos theta.x2) in
      Maths.[| [| F _A1 + (F 2.0 * z); F _A3 + z |]; [| F _A3 + z; F _A3 |] |]
      |> Maths.of_arrays
    in
    let c =
      let z = Maths.(F _A2 * sin theta.x2) in
      Maths.
        [| [| z * neg theta.x2_dot * ((F 2.0 * theta.x1_dot) + theta.x2_dot) |]
         ; [| z * sqr theta.x1_dot |]
        |]
      |> Maths.of_arrays
    in
    let h = Maths.of_arrays [| [| theta.x1_dot |]; [| theta.x2_dot |] |] in
    Linalg.linsolve m Maths.(torque - c - (b *@ h))


(* compute hand position from given angular state *)
let hand_of theta =
  let open Algodiff.D.Maths in
  let joint_x = F _L1 * cos theta.x1 in
  let joint_y = F _L1 * sin theta.x1 in
  let joint_x_dot = neg (F _L1 * theta.x1_dot * sin theta.x1) in
  let joint_y_dot = F _L1 * theta.x1_dot * cos theta.x1 in
  let z = theta.x1 + theta.x2 in
  let zdot = theta.x1_dot + theta.x2_dot in
  { x1 = joint_x + (F _L2 * cos z)
  ; x2 = joint_y + (F _L2 * sin z)
  ; x1_dot = joint_x_dot - (F _L2 * zdot * sin z)
  ; x2_dot = joint_y_dot + (F _L2 * zdot * cos z)
  }


(* given a time series of muscle activation, compute the hand trajectory *)
let theta_evolution ~dt torque =
  let n_bins = 1 + Algodiff.D.row_num torque in
  let theta = Array.make n_bins _THETA0 in
  let _theta = ref _THETA0 in
  for t = 1 to n_bins - 1 do
    let theta_dot_dot =
      theta_dot_dot !_theta Algodiff.D.Maths.(get_row torque (pred t))
    in
    let prev = !_theta in
    (_theta
       := Algodiff.D.Maths.
            { x1 = prev.x1 + (F dt * prev.x1_dot)
            ; x2 = prev.x2 + (F dt * prev.x2_dot)
            ; x1_dot = prev.x1_dot + (F dt * get_item theta_dot_dot 0 0)
            ; x2_dot = prev.x2_dot + (F dt * get_item theta_dot_dot 1 0)
            });
    theta.(t) <- !_theta
  done;
  theta


(* function that we expose *)
let theta_trajectory ~dt torque =
  let torque = Algodiff.D.pack_arr torque in
  theta_evolution ~dt torque


(* ----------------------------------------------------------------------------
   ---     TARGET STRAIGHT REACHES                                          ---
   ---------------------------------------------------------------------------- *)

let speed_profile =
  let tau_reach = 0.140 in
  let bell t =
    let t = t /. tau_reach in
    Maths.(sqr t *. exp (-.sqr t /. 2.))
  in
  let z = bell (tau_reach *. sqrt 2.) in
  fun t -> 1. /. z *. bell t


let _HAND0 = hand_of _THETA0

let straight_reach ~dt ~duration ~angle ~radius =
  let n_bins = Maths.(round (duration /. dt) |> int_of_float) in
  let ca = cos angle in
  let sa = sin angle in
  let hand_traj peak_speed =
    let hand = ref _HAND0 in
    Array.init n_bins (fun t ->
        let radial_speed = peak_speed *. speed_profile (dt *. float t) in
        let x1_dot = ca *. radial_speed in
        let x2_dot = sa *. radial_speed in
        let prev = !hand in
        let open Algodiff.D.Maths in
        hand
          := { x1 = prev.x1 + (F dt * prev.x1_dot)
             ; x2 = prev.x2 + (F dt * prev.x2_dot)
             ; x1_dot = F x1_dot
             ; x2_dot = F x2_dot
             };
        prev)
  in
  (* binary search to adjust peak speed to get to a desired radius *)
  let rec search_peak_speed a b =
    let middle = (a +. b) /. 2. in
    let hand_traj = hand_traj middle in
    let final = hand_traj.(pred (Array.length hand_traj)) in
    let achieved_radius =
      let dx1 = unpack final.x1 -. unpack _HAND0.x1
      and dx2 = unpack final.x2 -. unpack _HAND0.x2 in
      sqrt ((dx1 *. dx1) +. (dx2 *. dx2))
    in
    if abs_float (a -. b) < 0.0001
    then hand_traj, middle
    else if achieved_radius < radius
    then search_peak_speed middle b
    else search_peak_speed a middle
  in
  search_peak_speed 0. 3.


(* ----------------------------------------------------------------------------
   ---     OPTIMAL TORQUE                                                   ---
   ---------------------------------------------------------------------------- *)

(* get the optimal torque for a given target hand trajectory *)
let optimal_torque ~dt target_hands =
  let n_bins = Array.length target_hands in
  (* function that LBFGS will understand
     here prms and g are Array1.t *)
  let f _ prms =
    let _torque = prms.torques in
    let hands = Array.map hand_of (theta_evolution ~dt _torque) in
    let _, error =
      Array.fold_left
        (fun (i, accu) hi ->
          let thi = target_hands.(i) in
          let accu =
            Algodiff.D.Maths.(
              accu
              + sqr (hi.x1 - thi.x1)
              + sqr (hi.x1_dot - thi.x1_dot)
              + sqr (hi.x2 - thi.x2)
              + sqr (hi.x2_dot - thi.x2_dot))
          in
          succ i, accu)
        (0, F 0.)
        hands
    in
    let roughness =
      let pos = Algodiff.D.Maths.get_slice [ [ 1; pred (pred n_bins) ] ] _torque in
      let pre = Algodiff.D.Maths.get_slice [ [ 0; pred (pred (pred n_bins)) ] ] _torque in
      Algodiff.D.Maths.(l2norm_sqr' (pos - pre))
    in
    let init0 = Algodiff.D.Maths.(l2norm_sqr' (get_row _torque 0)) in
    Algodiff.D.Maths.(error + (F 10. * init0) + (F 2.0 * roughness))
  in
  (* random initial parameter vector *)
  let prms0 =
    { torques = Algodiff.D.(Maths.(Mat.(F 0.0001 * gaussian (pred n_bins) 2))) }
  in
  (* run BFGS *)
  let s = O.init ~prms0 () in
  let _ = O.min ~f s in
  let prms = O.prms s in
  Algodiff.D.unpack_arr prms.torques
