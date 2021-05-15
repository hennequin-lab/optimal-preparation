open Printf
open Lib
module T = Tspecs.Default

let () = Owl_stats_prng.init 0
let spontaneous, xstars, c = Simple.read_setup ~dynamic:false
let q = Simple.q_mov c
let qs = [| q |]
let xs, _ = Simple.naive ~xstars (module T)

let () =
  (* simulate the dynamics for each xstar *)
  Array.iteri
    (fun j xs ->
      Array.iteri
        (fun i x ->
          let xstar = xstars.(j).(i) in
          let path s = Simple.in_dir (sprintf "%s_naive_r%i_%i" s (succ j) (i + 1)) in
          Simple.save ~xstar ~qs (module T) ~path x)
        xs)
    xs
