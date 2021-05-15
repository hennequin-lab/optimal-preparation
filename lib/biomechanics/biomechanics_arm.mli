open Owl

type position
type angular

type 'a state =
  { x1 : Algodiff.D.t
  ; x1_dot : Algodiff.D.t
  ; x2 : Algodiff.D.t
  ; x2_dot : Algodiff.D.t
  }

(* takes a value in degrees, returns radians *)
val deg : float -> float
val unpack_sequence : 'a state array -> Mat.mat
val hand_of : angular state -> position state

(* angle: in radians;
   radius: in meters *)
val straight_reach
  :  dt:float
  -> duration:float
  -> angle:float
  -> radius:float
  -> position state array * float

(* computes the optimal torque for a given hand position trajectory;
   result is a time x 2 matrix *)
val optimal_torque : dt:float -> position state array -> Mat.mat
val theta_trajectory : dt:float -> Mat.mat -> angular state array
