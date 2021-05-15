open Owl
open Printf
open Lib

let () = Owl_stats_prng.init 0
let spontaneous, xstars, c = Simple.read_setup ~dynamic:false
let q = Simple.q_mov c

module T = Tspecs.Default

let xs, lqr_feedback = Simple.vanilla ~xstars ~q (module T)

let () =
  let qs = [| q |] in
  let yx, xy, _ = lqr_feedback |> Option.get in
  if not Cmdargs.(check "-reuse")
  then (
    Mat.save_txt yx ~out:Dir.(in_dir "vanilla_feedback_yx");
    Mat.save_txt xy ~out:Dir.(in_dir "vanilla_feedback_xy"));
  (* simulate the dynamics for each xstar *)
  Array.iteri
    (fun j xs ->
      Array.iteri
        (fun i (x, _, _) ->
          let xstar = xstars.(j).(i) in
          let path s =
            Simple.in_dir (sprintf "%s_vanilla_loop_r%i_%i" s (succ j) (succ i))
          in
          Simple.save ~qs ~xstar (module T) ~path (x, None, None);
          (* save the input u relative to the naive input *)
          Mat.save_txt Mat.(xy *@ yx *@ relu (x - xstar) |> transpose) ~out:(path "kx"))
        xs)
    xs
