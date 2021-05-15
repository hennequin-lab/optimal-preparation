module L = Lazy
open Owl

(* ----------------------------------------------------------------------------
   ---                   MISCELLANEOUS FUNCTIONS                            ---
   ---------------------------------------------------------------------------- *)

let save_bin filename m =
  let output = open_out filename in
  Marshal.to_channel output m [ Marshal.No_sharing ];
  close_out output


(* reads whatever was saved using [save_bin] *)
let read_bin filename =
  if not Sys.(file_exists filename)
  then failwith Printf.(sprintf "%s does not exist" filename)
  else (
    let input = open_in filename in
    let m = Marshal.from_channel input in
    close_in input;
    m)


let identity x = x
let delta i j = if i = j then 1. else 0.

let unit_trace q =
  let n, _ = Mat.shape q in
  let z = float n /. Mat.trace q in
  Mat.(z $* q)


let random_rotation x =
  let _, m = Mat.shape x in
  let rot, _, _ = Linalg.D.qr (Mat.gaussian m m) in
  Mat.(x *@ rot)


let random_h n desired_norm =
  let h = Mat.gaussian n 1 in
  Mat.(desired_norm /. l2norm' h $* h)


let all_but i m =
  let n = Array.length m in
  let ids =
    Array.append (Array.init i (fun k -> k)) (Array.init (n - i - 1) (fun k -> k + i + 1))
  in
  Array.to_list ids


let spectral_abscissa m = Linalg.D.eig m |> snd |> Dense.Matrix.Z.re |> Mat.max'

let eigenvalues m =
  let v = Linalg.D.eigvals m in
  let re = Dense.Matrix.Z.re v in
  let im = Dense.Matrix.Z.im v in
  Mat.(concat_horizontal (transpose re) (transpose im))


let maybe x f alt =
  match x with
  | Some x -> f x
  | None -> alt


let maybe_reuse reuse file x =
  if reuse
  then read_bin (Dir.in_dir file)
  else (
    let x = L.force x in
    save_bin (Dir.in_dir file) x;
    x)


let low_pass_filter ?init ~dt ~tau x =
  let _, n_bins = Mat.shape x in
  let ystore = Mat.(zeros (row_num x) n_bins) in
  let y =
    ref
      (match init with
      | Some v -> v
      | None -> Mat.zeros Mat.(row_num x) 1)
  in
  Mat.iteri_cols
    (fun t x ->
      Mat.(ystore.${[]; [ t ]} <- !y);
      y := Mat.((1. -. (dt /. tau) $* !y) + (dt /. tau $* x)))
    x;
  ystore
