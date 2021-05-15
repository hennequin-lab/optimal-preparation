(* --------------------------------------------------------------------------------
   --     GOOD DEFAULT PARAMETERS                                                --
   -------------------------------------------------------------------------------- *)

module Base = struct
  let dt = 2E-4
  let sampling_dt = 1E-3
  let tau = 150E-3
  let noise_prms = None
end

let baseline_rate = 1.
let spontaneous_std = 0.15
let mov_input = Signals.alpha_bump ~tau_rise:50E-3 ~tau_decay:0.5 ~amplitude:5.

(* energy costs *)
let r2_vanilla = 0.1
let r2_instant = 0.1
let r2_dynamic = 0.01
let taus = 10E-3, 10E-3

(* number of (phantom) muscles *)
let n_muscles_vanilla = 0
let n_muscles_instant = 0
let n_muscles_naive = 0
let n_muscles_dynamic = 0
let n_muscles_minimum_energy = 0

(* tau_z, tau_y *)
