open Owl
open Lib
open Miscellaneous
open Setup_base

(* do everything with the default set of parameters *)
open Defaults
open Defaults.Main

(* Gramians *)
module G = Gramians.Make (struct
  let a = a
  let m1_slice = m1_slice
  let gamma = None
end)

let n_obs = n
let path = Dir.in_dir
let n_mov = n_mov
let m1_slice = m1_slice
let n_muscles = n_muscles
let lambda_traj = 1.0 /. Int.to_float n_mov
let lambda_reg = 1.0 /. Int.to_float (List.length m1_slice * n_muscles)

let spontaneous =
  let get (a, _b, _c, _d) = a in
  if reuse
  then get (read_bin (path "setup.bin"))
  else (
    (* make up a new random one *)
    let v = Mat.gaussian n 1 in
    let v = Mat.(spontaneous_std /. std' v $* v) in
    Mat.(baseline_rate $+ v))


(* ----------- OPTIMIZATION OF XSTARS AND C ------------ *)
module O = Owl_opt_lbfgs.D.Make (Setup_base.Prms)
open Setup_base.Prms

let cost_tracker = Tracker.init Dir.(in_dir "cost_tracker")
let __w_rec = Algodiff.D.Arr w_rec
let __gamma = Algodiff.D.Arr gamma
let __spontaneous = Algodiff.D.Arr spontaneous
let __target_torques = time_array target_torques

let unpack =
  (* we'll assemble xstar using the top eigenvectors of the observability Gramian *)
  let top_obs = G.O.top n_obs in
  let z = Owl.Maths.(float n *. float n_mov *. sqr xstars_std) in
  (* unpack parameters *)
  fun prms ->
    let xstars_prms = prms.xstars
    and c = prms.c in
    let open Algodiff.D in
    let xstars = Maths.(Arr top_obs *@ xstars_prms) in
    let xstars = Maths.(sqrt (F z / l2norm_sqr' xstars) * xstars) in
    (* make sure C has all the xstars and spontaneous in its nullspace *)
    let xstars_motor =
      let z = Maths.(concatenate ~axis:1 [| xstars + __spontaneous; __spontaneous |]) in
      Maths.(__gamma *@ z)
    in
    let xstars_motor_t = Maths.(transpose xstars_motor) in
    let h = Linalg.linsolve Maths.(xstars_motor_t *@ xstars_motor) xstars_motor_t in
    let c = Maths.(c - (c *@ xstars_motor *@ h)) in
    Maths.(__spontaneous + xstars), c


let trajectory =
  let module DDiff =
    Dynamics_autodiff.Make (struct
      let dt = dt
      let sampling_dt = sampling_dt
      let tau = tau
      let spontaneous = spontaneous
      let mov_input = mov_input
    end)
  in
  let w_rec = Algodiff.D.pack_arr w_rec in
  fun xstars c ->
    let open Algodiff.D in
    let r = DDiff.run ~w_rec ~n_bins:it_max ~layers:None xstars in
    let t =
      let cgamma = Maths.(c *@ __gamma) in
      Array.map (fun r -> Maths.(cgamma *@ r)) r
    in
    r, t


let f =
  let rescaling = Mat.of_arrays [| [| 1. |]; [| 3. |] |] in
  fun _ prms ->
    let xstars, c = unpack prms in
    let open Algodiff.D in
    let _, torques = trajectory xstars c in
    let cost_move =
      let _, tmp =
        Array.fold_left
          (fun (t, accu) torque ->
            ( succ t
            , Maths.(
                accu + l2norm_sqr' (Arr rescaling * (torque - Arr __target_torques.(t))))
            ))
          (0, F 0.)
          torques
      in
      Maths.(F (lambda_traj *. sampling_dt) * tmp)
    in
    let reg = Maths.(F lambda_reg * l2norm_sqr' c) in
    Maths.(cost_move + reg)


let stop =
  let cost_move_old = ref 1E9 in
  let sa = Miscellaneous.(spectral_abscissa w_rec) in
  fun _ s ->
    let k = O.iter s in
    if k < 1
    then false
    else (
      let cost_value = O.prev_fv s in
      Tracker.update cost_tracker cost_value;
      let prms = O.prms s in
      let xstars, c = unpack prms in
      let c = AD.unpack_arr c in
      let cost_move = cost_value -. (lambda_reg *. Mat.(l2norm_sqr' c)) in
      if pred k mod 10 = 0
      then (
        let xstars = AD.unpack_arr xstars in
        Mat.save_txt c ~out:(path "c");
        Mat.save_txt spontaneous ~out:(path "spontaneous");
        Mat.save_txt xstars ~out:(path "xstars");
        Mat.save_txt Mat.(xstars - spontaneous) ~out:(path "dxstars");
        Mat.save_txt c ~out:(path "c");
        Mat.save_txt Mat.(c *@ gamma *@ xstars) ~out:(path "cxstars");
        Mat.save_txt Mat.(c *@ gamma *@ spontaneous) ~out:(path "csp");
        let rs, torques = trajectory (Arr xstars) (Arr c) in
        save_trajs_and_compute_hands ~path rs torques |> ignore;
        let xstars =
          Array.init n_radius (fun i ->
              Array.init n_angles (fun j -> Mat.col xstars (j + (n_angles * i))))
        in
        Tracker.flush cost_tracker;
        save_bin (path "setup.bin") (spontaneous, target_torques, xstars, c);
        save_bin (path "setup_prms.bin") prms);
      let pct_change = (!cost_move_old -. cost_move) /. !cost_move_old |> Float.abs in
      cost_move_old := cost_move;
      Stdio.printf
        "\riter %5i | cost = %.8f | sa %.3f | cost_move %f | pct_change %f %!"
        k
        cost_value
        sa
        cost_move
        pct_change;
      k >= max_iter || cost_move < 5E-4 || pct_change < 1E-4)


(* finish off with BFGS *)
let prms =
  (* xstars_prms : n_obs * n_mov
     c_prms : n_muscles * n_m1 *)
  let prms0 =
    if reuse
    then read_bin (path "setup_prms.bin")
    else
      { xstars = Algodiff.D.Mat.gaussian ~sigma:(0.1 /. sqrt (float n)) n_obs n_mov
      ; c =
          Algodiff.D.Mat.gaussian
            ~sigma:(0.1 /. sqrt (float n))
            n_muscles
            List.(length m1_slice)
      }
  in
  let s = O.init ~prms0 () in
  O.min ~stop ~f s |> ignore;
  O.prms s


let () =
  let xstars, c = unpack prms in
  let xstars = Algodiff.D.unpack_arr xstars in
  let c = Algodiff.D.unpack_arr c in
  Mat.save_txt c ~out:(path "c");
  Mat.save_txt spontaneous ~out:(path "spontaneous");
  Mat.save_txt xstars ~out:(path "xstars");
  let xstars =
    Array.init n_radius (fun i ->
        Array.init n_angles (fun j -> Mat.col xstars (j + (n_angles * i))))
  in
  save_bin (path "setup.bin") (spontaneous, target_torques, xstars, c)
