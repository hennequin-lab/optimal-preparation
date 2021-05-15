open Printf
open Owl
open Lib
open Biomechanics

(* invert the biomechanics of the arm + muscle model
   to get "muscle activations" for each of 8 straight reaches *)

let dt = 1E-3
let duration = 0.6
let n_radius = Cmdargs.(get_int "-n_radius" |> default 1)
let n_angles = Cmdargs.(get_int "-n_angles" |> default 8)
let reach_angles = Task.(reach_angles n_angles)
let radius = if n_radius = 1 then [| 0.2 |] else Task.(reach_radius n_radius)

let _ =
  Miscellaneous.save_bin Dir.(in_dir "reaches_radius.bin") radius;
  Miscellaneous.save_bin Dir.(in_dir "reaches_angles.bin") reach_angles;
  let open Arm in
  Array.mapi
    (fun jd radius ->
      let jd = succ jd in
      Array.mapi
        (fun id reach_angle ->
          print_newline ();
          let id = succ id in
          let hand_trajectory, _ =
            straight_reach ~dt ~duration ~angle:(deg reach_angle) ~radius
          in
          (* save the trajectory *)
          (hand_trajectory
          |> unpack_sequence
          |> fun m ->
          Mat.save_txt m ~out:(Dir.in_dir (sprintf "target_hand_traj_r%i_%i" jd id)));
          let torque = optimal_torque ~dt hand_trajectory in
          Miscellaneous.save_bin
            (Dir.in_dir (sprintf "target_torque_r%i_%i.bin" jd id))
            torque;
          Mat.save_txt torque ~out:(Dir.in_dir (sprintf "target_torque_r%i_%i" jd id));
          (* check if the torque obtained gives us the right trajectory *)
          let theta_trajectory = theta_trajectory ~dt torque in
          Mat.save_txt
            (unpack_sequence theta_trajectory)
            ~out:(Dir.in_dir (sprintf "actual_theta_traj_r%i_%i" jd id));
          let actual_hand_trajectory = Array.map hand_of theta_trajectory in
          Mat.save_txt
            (unpack_sequence actual_hand_trajectory)
            ~out:(Dir.in_dir (sprintf "actual_hand_traj_r%i_%i" jd id));
          torque)
        reach_angles)
    radius
  |> Miscellaneous.save_bin Dir.(in_dir "target_torques.bin")
