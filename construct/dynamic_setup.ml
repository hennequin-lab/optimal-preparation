open Printf
open Owl
module AD = Algodiff.D
open Lib
open Miscellaneous
open Biomechanics
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
let n_mov = n_mov
let m1_slice = m1_slice
let n_muscles = n_muscles
let lambda_traj = 1.0 /. Int.to_float n_mov
let lambda_reg = 1.0 /. Int.to_float (List.length m1_slice * n_muscles)

let spontaneous, c =
  let get (a, _b, _c, d) = a, d in
  if reuse
  then get (read_bin Dir.(in_dir "dynamic_setup.bin"))
  else get (read_bin Dir.(in_dir "setup.bin"))


let path s = Dir.(in_dir (sprintf "dynamic_%s" s))

(* ----------- OPTIMIZATION OF XSTARS AND C ------------ *)
module O = Owl_opt_lbfgs.D.Make (Setup_base.DynamicPrms)
open Setup_base.DynamicPrms

let __c = AD.Arr c
let __w_rec = AD.Arr w_rec
let __gamma = AD.Arr gamma
let __cgamma = AD.Maths.(__c *@ __gamma)
let __spontaneous = AD.Arr spontaneous
let __target_torques = time_array target_torques

let __top_obs =
  let top_obs = G.O.top n_obs in
  AD.Arr top_obs


let __cgamma_top_obs = AD.Maths.(__cgamma *@ __top_obs)

(* used to make sure that xstars stay in the null space of fixed c *)
let __projection =
  let open AD in
  let e = Mat.eye n_obs in
  let h =
    Linalg.linsolve
      Maths.(__cgamma_top_obs *@ transpose __cgamma_top_obs)
      __cgamma_top_obs
  in
  Maths.(e - (transpose __cgamma_top_obs *@ h))


let unpack =
  (* we'll assemble xstar using the top eigenvectors of the observability Gramian *)
  let z = Owl.Maths.(float n *. float n_mov *. sqr xstars_std) in
  (* unpack parameters *)
  fun prms ->
    let xstars_prms = prms.xstars in
    let open AD in
    let xstars = Maths.(__top_obs *@ __projection *@ xstars_prms) in
    let xstars = Maths.(sqrt (F z / l2norm_sqr' xstars) * xstars) in
    Maths.(__spontaneous + xstars)


let trajectory =
  let layers =
    let tau_z, tau_y = taus in
    let q = Simple.q_mov c in
    let _, (xz, zy, yx) =
      Miscellaneous.maybe_reuse
        reuse
        "dynamic_feedback.bin"
        (lazy
          (let feedback =
             Controller.dynamic ~a ~pmd_slice ~n_e ~tau ~taus ~r2:r2_dynamic ~q
           in
           let accuracy, decomposed =
             Controller.decompose ~n ~pmd_slice ~n_z:n ~reg:0.1 feedback
           in
           printf "decomposition accuracy = %f\n%!" accuracy;
           feedback, decomposed))
    in
    fun xstars ->
      let open AD in
      let z0 = Maths.(Arr zy *@ relu (Arr yx *@ relu xstars)) in
      let y0 = Maths.(Arr yx *@ relu xstars) in
      Some ((tau_z, xz, z0), (tau_y, zy, y0))
  in
  let module DDiff =
    Dynamics_autodiff.Make (struct
      let dt = dt
      let sampling_dt = sampling_dt
      let tau = tau
      let spontaneous = spontaneous
      let mov_input = mov_input
    end)
  in
  let w_rec = AD.pack_arr w_rec in
  fun xstars c ->
    let open AD in
    let r = DDiff.run ~w_rec ~n_bins:it_max ~layers:(layers xstars) xstars in
    let t =
      let cgamma = Maths.(c *@ __gamma) in
      Array.map (fun r -> Maths.(cgamma *@ r)) r
    in
    r, t


let f =
  let rescaling = Mat.of_arrays [| [| 1. |]; [| 3. |] |] in
  fun _ prms ->
    let xstars = unpack prms in
    let open AD in
    let _, torques = trajectory xstars __c in
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
    let reg = Maths.(F lambda_reg * l2norm_sqr' __c) in
    Maths.(cost_move + reg)


let stop =
  let cv = ref 1e9 in
  fun _ s ->
    let k = O.iter s in
    if k < 1
    then false
    else (
      let cost_value = O.prev_fv s in
      let pct_change = abs_float ((cost_value -. !cv) /. !cv) in
      cv := cost_value;
      printf "\riter %5i | cost = %.5f | pct change = %.5f%!" k cost_value pct_change;
      if pred k mod 10 = 0
      then (
        let prms = O.prms s in
        let xstars = unpack prms in
        let xstars = AD.unpack_arr xstars in
        Mat.save_txt c ~out:(path "c");
        Mat.save_txt spontaneous ~out:(path "spontaneous");
        Mat.save_txt xstars ~out:(path "xstars");
        Mat.save_txt Mat.(xstars - spontaneous) ~out:(path "dxstars");
        Mat.save_txt c ~out:(path "c");
        Mat.save_txt Mat.(c *@ gamma *@ xstars) ~out:(path "cxstars");
        Mat.save_txt Mat.(c *@ gamma *@ spontaneous) ~out:(path "csp");
        let rs, torques = trajectory (Arr xstars) (Arr c) in
        let rs = Array.map AD.unpack_arr rs |> mov_array in
        let torques = Array.map AD.unpack_arr torques |> mov_array in
        let thetas =
          Array.map
            (Array.map (fun torque ->
                 Arm.theta_trajectory ~dt:sampling_dt (Mat.transpose torque)))
            torques
        in
        let hands = Array.map (Array.map (Array.map Arm.hand_of)) thetas in
        Array.iteri
          (fun j thetas ->
            Array.iteri
              (fun k theta ->
                Arm.unpack_sequence theta
                |> Mat.(save_txt ~out:(path (sprintf "theta_r%i_%i" (succ j) (succ k)))))
              thetas)
          thetas;
        Array.iteri
          (fun j hands ->
            Array.iteri
              (fun k hand ->
                Arm.unpack_sequence hand
                |> Mat.(save_txt ~out:(path (sprintf "hand_r%i_%i" (succ j) (succ k)))))
              hands)
          hands;
        Array.iteri
          (fun j rs ->
            Array.iteri
              (fun k m ->
                Mat.transpose m
                |> Mat.(save_txt ~out:(path (sprintf "x_r%i_%i" (succ j) (succ k)))))
              rs)
          rs;
        Array.iteri
          (fun j torques ->
            Array.iteri
              (fun k m ->
                Mat.(transpose m)
                |> Mat.(save_txt ~out:(path (sprintf "torque_r%i_%i" (succ j) (succ k)))))
              torques)
          torques;
        let xstars =
          Array.init n_radius (fun i ->
              Array.init n_angles (fun j -> Mat.col xstars (j + (n_angles * i))))
        in
        save_bin Dir.(in_dir "dynamic_setup.bin") (spontaneous, target_torques, xstars, c);
        save_bin Dir.(in_dir "dynamic_setup_prms.bin") prms);
      k >= max_iter || pct_change > 1e-3)


(* finish off with BFGS *)
let prms =
  (* xstars_prms : n_obs * n_mov
     c_prms : n_muscles * n_m1 *)
  let prms0 =
    if reuse
    then read_bin Dir.(in_dir "dynamic_setup_prms.bin")
    else (
      let prms = read_bin Dir.(in_dir "setup_prms.bin") in
      let xstars_prms = Setup_base.Prms.(prms.xstars) in
      { xstars = AD.Linalg.(linsolve __projection xstars_prms) })
  in
  let s = O.init ~prms0 () in
  O.min ~stop ~f s |> ignore;
  O.prms s


let () =
  let xstars = unpack prms in
  let xstars = AD.unpack_arr xstars in
  Mat.save_txt c ~out:(path "c");
  Mat.save_txt spontaneous ~out:(path "spontaneous");
  Mat.save_txt xstars ~out:(path "xstars");
  let xstars =
    Array.init n_radius (fun i ->
        Array.init n_angles (fun j -> Mat.col xstars (j + (n_angles * i))))
  in
  save_bin Dir.(in_dir "dynamic_setup.bin") (spontaneous, target_torques, xstars, c)
