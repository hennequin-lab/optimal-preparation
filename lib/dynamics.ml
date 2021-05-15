open Printf
open Owl
module P = Stdlib

module type PT = sig
  val n : int
  val n_e : int
  val n_i : int
  val w_rec : Mat.mat
  val spontaneous : Mat.mat
  val dt : float
  val tau : float
  val nl : Mat.mat -> Mat.mat

  (* time constant of cortical dynamics *)
  val noise_prms : (float * Mat.mat * (float -> float)) option

  (* optional input noise parameters *)
end

module Make (P : PT) = struct
  include P

  let a = Mat.(w_rec - eye n)
  let h = Mat.(spontaneous - (w_rec *@ nl spontaneous))

  let prepare_noise ?eta0 n_trials =
    match noise_prms with
    | None -> None
    | Some (tau_eta, ell, noise_gain) ->
      let m = Mat.col_num ell in
      let eta =
        match eta0 with
        | Some eta0 -> ref eta0
        | None -> ref Mat.(ell *@ gaussian m n_trials)
      in
      let decay = 1. -. (dt /. tau_eta) in
      let factor =
        let sqrtdt = sqrt dt in
        sqrtdt /. tau_eta
      in
      Some
        (fun time ->
          (eta := Mat.((decay $* !eta) + (factor $* ell *@ gaussian m n_trials)));
          Mat.(noise_gain time $* !eta))


  (* simulation with instantaneous loops *)
  let simulate_instantaneous_feedback ?eta0 ~sampling_dt ~duration ~mov_input ~x0 =
    let x0, n_trials =
      match x0 with
      | Some v -> Mat.copy v, snd (Mat.shape v)
      | None -> Mat.zeros n 1, 1
    in
    let n_bins = Maths.round (duration /. dt) |> int_of_float in
    let sample_every = Maths.round (sampling_dt /. dt) |> int_of_float in
    let noise = prepare_noise ?eta0 n_trials in
    fun ?(verbose = false) ?perturbation ?(nl = nl) ~loop xstar ->
      let (yx, xy), loop_gain = loop in
      let feedback = Mat.(xy *@ yx) in
      let specific_input = Mat.(xstar - h - ((w_rec + feedback) *@ nl xstar)) in
      let xstore = Arr.empty [| n_bins / sample_every; n; n_trials |] in
      let x = ref x0 in
      for t = 0 to n_bins - 1 do
        let time = dt *. float t in
        if t mod 10 = 0
        then (
          if verbose then printf "\r[dynamics] %.3f / %.3f%!" time duration;
          Gc.minor ());
        let r = nl !x in
        let input =
          Mat.(
            h
            +$ mov_input time
            + (loop_gain time $* specific_input)
            + (w_rec *@ r)
            + (loop_gain time $* feedback *@ r))
        in
        let input =
          match noise with
          | Some draw -> Mat.(input + draw time)
          | None -> input
        in
        let input =
          match perturbation with
          | Some p -> Mat.(input + p time)
          | None -> input
        in
        (x := Mat.((P.(1. -. (dt /. tau)) $* !x) + (P.(dt /. tau) $* input)));
        if t mod sample_every = 0
        then Arr.(copy_ r ~out:(slice_left xstore [| Stdlib.(t / sample_every) |]))
      done;
      if verbose then print_newline ();
      let result = Arr.transpose ~axis:[| 2; 1; 0 |] xstore (* trials; neurons; time *) in
      (* squeeze the result if n_trials = 1 *)
      if n_trials = 1 then Arr.squeeze result else result


  (* simulation with dynamic loops *)
  let simulate_dynamic_feedback ?eta0 ~sampling_dt ~duration ~mov_input ~x0 =
    let x0, n_trials =
      match x0 with
      | Some v -> Mat.copy v, snd (Mat.shape v)
      | None -> Mat.zeros n 1, 1
    in
    let n_bins = Maths.round (duration /. dt) |> int_of_float in
    let sample_every = Maths.round (sampling_dt /. dt) |> int_of_float in
    let noise = prepare_noise ?eta0 n_trials in
    fun ?(verbose = false) ?perturbation ?(nl = nl) ~taus ~loop xstar ->
      let (xz, zy, yx), loop_gain, input_gain = loop in
      let feedback = Mat.(xz *@ zy *@ yx) in
      let specific_input = Mat.(xstar - h - ((w_rec + feedback) *@ nl xstar)) in
      let xstore = Arr.empty [| n_bins / sample_every; n; n_trials |] in
      let ystore = Arr.empty [| n_bins / sample_every; Mat.row_num yx; n_trials |] in
      let zstore = Arr.empty [| n_bins / sample_every; Mat.row_num zy; n_trials |] in
      let tau_z, tau_y = taus in
      let x = ref x0 in
      (* initialised in their steady-state given cortex *)
      let y = ref Mat.(zeros (row_num yx) n_trials) in
      let z = ref Mat.(zeros (row_num zy) n_trials) in
      (* run the dynamics forward in time *)
      for t = 0 to n_bins - 1 do
        let time = dt *. float t in
        if t mod 10 = 0
        then (
          Gc.minor ();
          if verbose then printf "\r[dynamics] %.3f / %.3f%!" time duration);
        let rx = nl !x
        and ry = !y
        and rz = nl !z in
        (* update cortex *)
        let input =
          Mat.(
            h
            +$ mov_input time
            + (w_rec *@ rx)
            + (xz *@ rz)
            + (input_gain time $* specific_input))
        in
        let input =
          match noise with
          | Some draw -> Mat.(input + draw time)
          | None -> input
        in
        let input =
          match perturbation with
          | Some p -> Mat.(input + p time)
          | None -> input
        in
        (x := Mat.((P.(1. -. (dt /. tau)) $* !x) + (P.(dt /. tau) $* input)));
        (* update thal + inter *)
        (z := Mat.((P.(1. -. (dt /. tau_z)) $* !z) + (P.(dt /. tau_z) $* zy *@ ry)));
        (y
           := Mat.(
                (P.(1. -. (dt /. tau_y)) $* !y)
                + (P.(loop_gain time *. dt /. tau_y) $* yx *@ rx)));
        (* store activity in buffer *)
        if t mod sample_every = 0
        then (
          let s = [| t / sample_every |] in
          Arr.(copy_ rx ~out:(slice_left xstore s));
          Arr.(copy_ ry ~out:(slice_left ystore s));
          Arr.(copy_ rz ~out:(slice_left zstore s)))
      done;
      if verbose then print_newline ();
      let reformat v =
        v
        |> Arr.transpose ~axis:[| 2; 1; 0 |]
        |> fun v -> if n_trials = 1 then Arr.squeeze v else v
      in
      reformat xstore, reformat ystore, reformat zstore
end
