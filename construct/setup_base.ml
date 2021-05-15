open Base
open Owl
open Lib
open Biomechanics
open Miscellaneous
open Defaults_base.Base
module AD = Algodiff.D

module Prms = struct
  type 'a t =
    { xstars : 'a
    ; c : 'a
    }
  [@@deriving prms]
end

module DynamicPrms = struct
  type 'a t = { xstars : 'a } [@@deriving prms]
end

let in_target_dir =
  let dir = Cmdargs.(get_string "-target_dir" |> default Dir.dir) in
  Printf.sprintf "%s/%s" dir


let max_iter = Cmdargs.(get_int "-max_iter" |> default 500000)
let reuse = Cmdargs.check "-reuse"

let target_torques =
  try read_bin (in_target_dir "target_torques.bin") with
  | _ ->
    Array.init 1 ~f:(fun i ->
        Array.init 8 ~f:(fun j ->
            Mat.load_txt
              (in_target_dir
                 Printf.(sprintf "target_torque_r%i_%i" Int.(succ i) Int.(succ j)))))


let n_radius = Array.length target_torques
let n_angles = Array.length target_torques.(0)
let n_mov = n_radius * n_angles
let n_muscles = snd (Mat.shape target_torques.(0).(0))
let xstars_std = 0.2
let n_bins = fst (Mat.shape target_torques.(0).(0))
let t_max = sampling_dt *. Int.to_float n_bins
let it_max = Float.to_int (t_max /. dt)
let n_bins = Float.to_int (t_max /. sampling_dt)
let sample_every = Float.to_int (sampling_dt /. dt)

(* reformat EMGs for the purpose of comparing with the data generated
   by the network below *)
let time_array x =
  Array.init n_bins ~f:(fun t ->
      Mat.concatenate
        ~axis:1
        (Array.map
           ~f:(fun y ->
             Mat.concatenate
               ~axis:1
               (Array.map ~f:(fun m -> Mat.transpose (Mat.row m t)) y))
           x))


let mov_array x =
  Array.init n_radius ~f:(fun y ->
      Array.init n_angles ~f:(fun m ->
          Mat.concatenate
            ~axis:1
            (Array.map ~f:(fun t -> Mat.col t (m + (y * n_angles))) x)))


let lambda_traj = 1.0 /. Int.to_float n_mov
let lambda_reg = 1.0 /. Int.to_float (List.length Defaults.Main.m1_slice * n_muscles)

module Tracker = struct
  type t =
    { history : float list ref
    ; buffer : float list ref
    ; file : string
    }

  let init file =
    (* clear file history *)
    Stdio.Out_channel.with_file file ~f:(fun _ -> ());
    { history = ref []; buffer = ref []; file }


  let update { history; buffer; file = _ } c =
    history := c :: !history;
    buffer := c :: !buffer


  let flush { history = _; buffer; file } =
    Stdio.Out_channel.with_file ~append:true file ~f:(fun out ->
        Stdio.Out_channel.output_lines out List.(rev_map ~f:Float.to_string !buffer));
    buffer := []
end

let save_trajs_and_compute_hands ~path rs torques =
  let rs = Array.map ~f:AD.unpack_arr rs |> mov_array in
  let torques = Array.map ~f:AD.unpack_arr torques |> mov_array in
  let thetas =
    Array.map
      ~f:
        (Array.map ~f:(fun torque ->
             Arm.theta_trajectory ~dt:sampling_dt (Mat.transpose torque)))
      torques
  in
  let hands =
    Array.map
      ~f:
        (Array.map ~f:(fun theta -> Array.map ~f:Arm.hand_of theta |> Arm.unpack_sequence))
      thetas
  in
  Array.iteri
    ~f:(fun j thetas ->
      Array.iteri
        ~f:(fun k theta ->
          Arm.unpack_sequence theta
          |> Mat.(
               save_txt
                 ~out:(path (Printf.sprintf "theta_r%i_%i" (Int.succ j) (Int.succ k)))))
        thetas)
    thetas;
  Array.iteri
    ~f:(fun j hands ->
      Array.iteri
        ~f:(fun k hand ->
          hand
          |> Mat.(
               save_txt
                 ~out:(path (Printf.sprintf "hand_r%i_%i" (Int.succ j) (Int.succ k)))))
        hands)
    hands;
  Array.iteri
    ~f:(fun j rs ->
      Array.iteri
        ~f:(fun k m ->
          Mat.transpose m
          |> Mat.(
               save_txt ~out:(path (Printf.sprintf "x_r%i_%i" (Int.succ j) (Int.succ k)))))
        rs)
    rs;
  Array.iteri
    ~f:(fun j torques ->
      Array.iteri
        ~f:(fun k m ->
          Mat.(transpose m)
          |> Mat.(
               save_txt
                 ~out:(path (Printf.sprintf "torque_r%i_%i" (Int.succ j) (Int.succ k)))))
        torques)
    torques
