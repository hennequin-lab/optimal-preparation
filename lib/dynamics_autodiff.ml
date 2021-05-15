(* dynamics of the movement phase, written under Algodiff
   to be used in the optimization of C and xstars *)

open Owl
module AD = Algodiff.D
module P = Stdlib

module type PT = sig
  val dt : float
  val sampling_dt : float
  val tau : float
  val spontaneous : Mat.mat
  val mov_input : float -> float
end

module Make (P : PT) = struct
  open P

  let sample_every = Maths.round (sampling_dt /. dt) |> int_of_float
  let spontaneous = AD.pack_arr spontaneous

  let run ?(nl = AD.Maths.relu) ~w_rec ~n_bins ~layers xstars =
    let open AD in
    let h = Maths.(spontaneous - (w_rec *@ nl spontaneous)) in
    let rec iter t x layers accu =
      if t = n_bins
      then accu |> List.rev |> Array.of_list
      else (
        let time = Stdlib.(dt *. float t) in
        (* update x *)
        let r = nl x in
        let accu = if t mod sample_every = 0 then r :: accu else accu in
        let input = Maths.(F (mov_input time) + h + (w_rec *@ r) - x) in
        (* perhaps add the input left over from (still decaying) previous layer z *)
        let input =
          match layers with
          | Some ((_, xz, z), _) -> Maths.(input + (Arr xz *@ nl z))
          | None -> input
        in
        let new_x = Maths.(x + (F Stdlib.(dt /. tau) * input)) in
        (* maybe udpate z *)
        let new_layers =
          match layers with
          | Some ((tau_z, xz, z), (tau_y, zy, y)) ->
            let y = Maths.(F Stdlib.(1. -. (dt /. tau_y)) * y) in
            let z =
              Maths.(
                (F Stdlib.(1. -. (dt /. tau_z)) * z) + (F (dt /. tau_z) * (Arr zy *@ y)))
            in
            Some ((tau_z, xz, z), (tau_y, zy, y))
          | None -> None
        in
        iter (succ t) new_x new_layers accu)
    in
    iter 0 xstars layers []
end
