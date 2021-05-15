open Owl

(* ----------------------------------------------------------------------------
   ---                   MISCELLANEOUS FUNCTIONS                            ---
   ---------------------------------------------------------------------------- *)

val identity : 'a -> 'a
val delta : int -> int -> float
val unit_trace : Mat.mat -> Mat.mat
val random_rotation : Mat.mat -> Mat.mat
val random_h : int -> float -> Mat.mat
val all_but : int -> 'a array -> int list
val spectral_abscissa : Mat.mat -> float
val eigenvalues : Mat.mat -> Mat.mat
val maybe : 'a option -> ('a -> 'b) -> 'b -> 'b
val maybe_reuse : bool -> string -> 'a lazy_t -> 'a
val low_pass_filter : ?init:Mat.mat -> dt:float -> tau:float -> Mat.mat -> Mat.mat
val save_bin : string -> 'a -> unit

(* reads whatever was saved using [save_bin] *)
val read_bin : string -> 'a
