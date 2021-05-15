open Owl
open Printf
open Lib

let save_hist w filename =
  let h = Stats.(histogram (`N 15) Mat.(to_array w) |> normalise_density) in
  let bins = Mat.of_array h.bins (-1) 1 in
  let bins =
    Mat.(0.5 $* get_slice [ [ 0; -2 ]; [] ] bins + get_slice [ [ 1; -1 ]; [] ] bins)
  in
  let d =
    match h.density with
    | Some d -> Mat.of_array d (-1) 1
    | None -> assert false
  in
  Mat.(
    save_txt
      (concatenate ~axis:1 [| bins; d; cumsum ~axis:0 d |])
      ~out:Dir.(in_dir filename))


module T = Tspecs.Default

let spontaneous, xstars, c = Simple.read_setup ~dynamic:true
let q = Simple.q_mov c
let qs = [| q |]
let activities, decomposed = Simple.dynamic ~xstars (module T)

let () =
  (* save the loop matrices *)
  let xz, zy, yx = decomposed |> Option.get in
  let yx = Option.get yx in
  if not Cmdargs.(check "-resue")
  then (
    Mat.save_txt xz ~out:Dir.(in_dir "dynamic_feedback_xz");
    Mat.save_txt zy ~out:Dir.(in_dir "dynamic_feedback_zy");
    Mat.save_txt yx ~out:Dir.(in_dir "dynamic_feedback_yx");
    (* save weight distributions *)
    save_hist xz "dynamic_feedback_xz_hist";
    save_hist zy "dynamic_feedback_zy_hist";
    save_hist
      Mat.(get_slice [ []; [ 0; pred Defaults.Main.n_e ] ] yx)
      "dynamic_feedback_yx_hist";
    save_hist Mat.(abs xz) "dynamic_feedback_xz_abs_hist";
    save_hist Mat.(abs zy) "dynamic_feedback_zy_abs_hist";
    save_hist
      Mat.(get_slice [ []; [ 0; pred Defaults.Main.n_e ] ] (abs yx))
      "dynamic_feedback_yx_abs_hist");
  (* simulate the dynamics for each xstar *)
  Array.iteri
    (fun j activities ->
      Array.iteri
        (fun i x ->
          let xstar = xstars.(j).(i) in
          let path s =
            Simple.in_dir (sprintf "%s_dynamic_loop_r%i_%i" s (succ j) (i + 1))
          in
          Simple.save ~qs ~xstar (module T) ~path x)
        activities)
    activities
